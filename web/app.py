#!/usr/bin/env python3
"""
PairTCR Web Application
A web interface for running the PairTCR pipeline
"""

import os
import sys
import time
import uuid
import shutil
import zipfile
import threading
import json
from datetime import datetime
from pathlib import Path
from flask import Flask, render_template, request, jsonify, send_file, url_for
from flask_socketio import SocketIO, emit, join_room, leave_room
import subprocess
import logging
from queue import Queue
import tempfile

# Add parent directory to Python path to import our scripts
web_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(web_dir)
sys.path.insert(0, parent_dir)

# Import our pipeline runner
from web.pipeline_runner import WebPipelineRunner

app = Flask(__name__)
app.config['SECRET_KEY'] = 'pairtcr-web-secret-key'
socketio = SocketIO(app, cors_allowed_origins="*")

# Global variables for job management
jobs = {}
job_queue = Queue()

class JobManager:
    def __init__(self):
        self.jobs = {}
        self.lock = threading.Lock()
    
    def create_job(self, job_id, params):
        with self.lock:
            self.jobs[job_id] = {
                'id': job_id,
                'status': 'pending',
                'created_at': datetime.now().isoformat(),
                'started_at': None,
                'completed_at': None,
                'params': params,
                'progress': 0,
                'current_step': '',
                'log_messages': [],
                'error_message': None,
                'result_path': None
            }
        return self.jobs[job_id]
    
    def update_job(self, job_id, **kwargs):
        with self.lock:
            if job_id in self.jobs:
                self.jobs[job_id].update(kwargs)
                # Emit update to websocket
                socketio.emit('job_update', self.jobs[job_id], room=job_id)
    
    def get_job(self, job_id):
        with self.lock:
            return self.jobs.get(job_id)
    
    def add_log_message(self, job_id, message):
        with self.lock:
            if job_id in self.jobs:
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                log_entry = f"[{timestamp}] {message}"
                self.jobs[job_id]['log_messages'].append(log_entry)
                # Emit log update to websocket
                socketio.emit('log_update', {
                    'job_id': job_id,
                    'message': log_entry
                }, room=job_id)

job_manager = JobManager()

def run_pipeline_job(job_id, params):
    """Run the PairTCR pipeline in a separate thread"""
    try:
        job_manager.update_job(job_id, 
                             status='running', 
                             started_at=datetime.now().isoformat(),
                             current_step='Initializing')
        
        job_manager.add_log_message(job_id, "Starting PairTCR pipeline...")
        job_manager.add_log_message(job_id, f"Parameters: {params}")
        
        # Determine output directory
        if params.get('output_dir') and params['output_dir'].strip():
            # User specified output directory
            output_root = params['output_dir'].strip()
            job_manager.add_log_message(job_id, f"Using user-specified output directory: {output_root}")
        else:
            # Default: same level as input directory with prefix name
            input_dir = params['input_dir']
            prefix = params.get('prefix', 'TCR_TSO_18')
            parent_dir_of_input = os.path.dirname(input_dir.rstrip('/'))
            output_root = os.path.join(parent_dir_of_input, prefix)
            job_manager.add_log_message(job_id, f"Using default output directory: {output_root}")
        
        # Create output directory
        os.makedirs(output_root, exist_ok=True)
        job_manager.add_log_message(job_id, f"Created output directory: {output_root}")
        
        # Create pipeline runner with web-specific logging
        pipeline = WebPipelineRunner(
            input_dir=params['input_dir'],
            output_root=output_root,
            prefix=params.get('prefix', 'TCR_TSO_18'),
            read_limit=params.get('read_limit', 100000),
            threads=params.get('threads', 4),
            mixcr_jar=params.get('mixcr_jar', os.path.join(parent_dir, 'scripts', 'mixcr.jar')),
            force_restart=params.get('force_restart', False),
            use_c_version=params.get('use_c_version', True),
            job_id=job_id,
            job_manager=job_manager
        )
        
        # Run the pipeline
        job_manager.update_job(job_id, current_step='Running pipeline', progress=10)
        pipeline.run_pipeline()
        
        # Create downloadable zip file
        job_manager.update_job(job_id, current_step='Creating download package', progress=90)
        job_manager.add_log_message(job_id, "Creating downloadable package...")
        
        # Create web results directory for ZIP files
        web_results_dir = os.path.join(web_dir, "results")
        os.makedirs(web_results_dir, exist_ok=True)
        zip_path = os.path.join(web_results_dir, f"{job_id}_results.zip")
        create_results_zip(output_root, zip_path)
        
        # Job completed successfully
        job_manager.update_job(job_id,
                             status='completed',
                             completed_at=datetime.now().isoformat(),
                             current_step='Completed',
                             progress=100,
                             result_path=zip_path)
        
        job_manager.add_log_message(job_id, "Pipeline completed successfully!")
        job_manager.add_log_message(job_id, f"Results available for download: {zip_path}")
        
    except Exception as e:
        error_msg = f"Pipeline failed: {str(e)}"
        job_manager.update_job(job_id,
                             status='failed',
                             completed_at=datetime.now().isoformat(),
                             current_step='Failed',
                             error_message=error_msg)
        job_manager.add_log_message(job_id, f"ERROR: {error_msg}")

def create_results_zip(source_dir, zip_path):
    """Create a zip file containing all results"""
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, source_dir)
                zipf.write(file_path, arcname)

@app.route('/')
def index():
    """Main page"""
    return render_template('index.html')

@app.route('/api/start_job', methods=['POST'])
def start_job():
    """Start a new pipeline job"""
    try:
        params = request.json
        
        # Validate required parameters
        if not params.get('input_dir'):
            return jsonify({'error': 'Input directory is required'}), 400
        
        if not os.path.exists(params['input_dir']):
            return jsonify({'error': f"Input directory '{params['input_dir']}' does not exist"}), 400
        
        # Generate unique job ID
        job_id = str(uuid.uuid4())
        
        # Create job
        job = job_manager.create_job(job_id, params)
        
        # Start pipeline in background thread
        thread = threading.Thread(target=run_pipeline_job, args=(job_id, params))
        thread.daemon = True
        thread.start()
        
        return jsonify({'job_id': job_id, 'status': 'started'})
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/job/<job_id>')
def get_job_status(job_id):
    """Get job status"""
    job = job_manager.get_job(job_id)
    if not job:
        return jsonify({'error': 'Job not found'}), 404
    return jsonify(job)

@app.route('/api/download/<job_id>')
def download_results(job_id):
    """Download job results"""
    job = job_manager.get_job(job_id)
    if not job:
        return jsonify({'error': 'Job not found'}), 404
    
    if job['status'] != 'completed':
        return jsonify({'error': 'Job not completed yet'}), 400
    
    if not job.get('result_path') or not os.path.exists(job['result_path']):
        return jsonify({'error': 'Results file not found'}), 404
    
    return send_file(job['result_path'], 
                     as_attachment=True, 
                     download_name=f"pairtcr_results_{job_id}.zip")

@app.route('/job/<job_id>')
def job_page(job_id):
    """Job monitoring page"""
    job = job_manager.get_job(job_id)
    if not job:
        return "Job not found", 404
    return render_template('job.html', job=job)

@socketio.on('join_job')
def on_join_job(data):
    """Join a job room for real-time updates"""
    job_id = data.get('job_id')
    if job_id:
        join_room(job_id)
        emit('joined', {'job_id': job_id})

@socketio.on('leave_job')
def on_leave_job(data):
    """Leave a job room"""
    job_id = data.get('job_id')
    if job_id:
        leave_room(job_id)

if __name__ == '__main__':
    # Create necessary directories
    os.makedirs(os.path.join(web_dir, 'results'), exist_ok=True)
    os.makedirs(os.path.join(web_dir, 'static'), exist_ok=True)
    
    # Run the Flask-SocketIO app
    print("Starting PairTCR Web Application...")
    print("Access the application at: http://localhost:5000")
    socketio.run(app, host='0.0.0.0', port=5000, debug=False) 