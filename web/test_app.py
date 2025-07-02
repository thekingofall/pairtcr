#!/usr/bin/env python3
"""
Simple test script for PairTCR Web Application
"""

import unittest
import tempfile
import os
import sys
import json
from pathlib import Path

# Add parent directory to path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

try:
    from app import app, job_manager
    FLASK_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import Flask app: {e}")
    FLASK_AVAILABLE = False

class TestPairTCRWebApp(unittest.TestCase):
    """Test cases for PairTCR Web Application"""
    
    def setUp(self):
        """Set up test fixtures"""
        if not FLASK_AVAILABLE:
            self.skipTest("Flask app not available")
        
        self.app = app
        self.app.config['TESTING'] = True
        self.client = self.app.test_client()
        
        # Create temporary directory for testing
        self.test_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures"""
        import shutil
        if hasattr(self, 'test_dir') and os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
    
    def test_index_page(self):
        """Test that index page loads"""
        response = self.client.get('/')
        self.assertEqual(response.status_code, 200)
        self.assertIn(b'PairTCR Pipeline', response.data)
    
    def test_start_job_missing_params(self):
        """Test starting job with missing parameters"""
        response = self.client.post('/api/start_job', 
                                  data=json.dumps({}),
                                  content_type='application/json')
        self.assertEqual(response.status_code, 400)
        data = json.loads(response.data)
        self.assertIn('error', data)
    
    def test_start_job_nonexistent_dir(self):
        """Test starting job with non-existent directory"""
        params = {
            'input_dir': '/nonexistent/directory',
            'prefix': 'test',
            'read_limit': 1000,
            'threads': 2
        }
        response = self.client.post('/api/start_job',
                                  data=json.dumps(params),
                                  content_type='application/json')
        self.assertEqual(response.status_code, 400)
        data = json.loads(response.data)
        self.assertIn('does not exist', data['error'])
    
    def test_job_manager(self):
        """Test job manager functionality"""
        # Create a test job
        job_id = 'test-job-123'
        params = {'input_dir': '/test', 'prefix': 'test'}
        
        job = job_manager.create_job(job_id, params)
        self.assertEqual(job['id'], job_id)
        self.assertEqual(job['status'], 'pending')
        
        # Update job
        job_manager.update_job(job_id, status='running', progress=50)
        updated_job = job_manager.get_job(job_id)
        self.assertEqual(updated_job['status'], 'running')
        self.assertEqual(updated_job['progress'], 50)
        
        # Add log message
        job_manager.add_log_message(job_id, 'Test message')
        job_with_logs = job_manager.get_job(job_id)
        self.assertEqual(len(job_with_logs['log_messages']), 1)
        self.assertIn('Test message', job_with_logs['log_messages'][0])

def test_imports():
    """Test that all required modules can be imported"""
    print("Testing imports...")
    
    required_modules = [
        'flask',
        'flask_socketio', 
        'pandas',
        'numpy',
        'tqdm'
    ]
    
    missing_modules = []
    for module in required_modules:
        try:
            __import__(module)
            print(f"✓ {module}")
        except ImportError:
            missing_modules.append(module)
            print(f"✗ {module}")
    
    if missing_modules:
        print(f"\nMissing modules: {missing_modules}")
        print("Install with: pip install -r requirements.txt")
        return False
    else:
        print("\n✓ All required modules available")
        return True

def test_file_structure():
    """Test that required files exist"""
    print("\nTesting file structure...")
    
    web_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(web_dir)
    
    required_files = [
        os.path.join(web_dir, 'app.py'),
        os.path.join(web_dir, 'pipeline_runner.py'),
        os.path.join(web_dir, 'requirements.txt'),
        os.path.join(web_dir, 'templates', 'index.html'),
        os.path.join(web_dir, 'templates', 'job.html'),
        os.path.join(parent_dir, 'scripts', 'pipeline_runner_base.py'),
    ]
    
    missing_files = []
    for file_path in required_files:
        if os.path.exists(file_path):
            print(f"✓ {os.path.relpath(file_path, web_dir)}")
        else:
            missing_files.append(file_path)
            print(f"✗ {os.path.relpath(file_path, web_dir)}")
    
    if missing_files:
        print(f"\nMissing files: {len(missing_files)}")
        return False
    else:
        print("\n✓ All required files found")
        return True

def main():
    """Run tests"""
    print("PairTCR Web Application Test Suite")
    print("=" * 50)
    
    # Test imports
    imports_ok = test_imports()
    
    # Test file structure  
    files_ok = test_file_structure()
    
    if not (imports_ok and files_ok):
        print("\n❌ Pre-requisites not met. Please fix the issues above.")
        return False
    
    # Run unit tests if Flask is available
    if FLASK_AVAILABLE:
        print("\nRunning unit tests...")
        unittest.main(verbosity=2, exit=False)
    else:
        print("\n⚠️ Flask app not available, skipping unit tests")
    
    print("\n✓ Test suite completed")
    return True

if __name__ == '__main__':
    main() 