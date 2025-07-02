#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PairTCR Web Server Startup Script
Checks dependencies and starts the web application
"""

import sys
import os
import subprocess
import argparse

def check_python_version():
    """Check if Python version is supported"""
    if sys.version_info < (3, 6):
        print("Error: Python 3.6 or higher is required")
        print("Current version: {}".format(sys.version))
        return False
    return True

def check_dependencies():
    """Check if required packages are installed"""
    required_packages = [
        'flask',
        'flask_socketio',
        'pandas',
        'tqdm',
        'numpy'
    ]
    
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        print("Missing required packages:")
        for package in missing_packages:
            print("  - {}".format(package))
        print()
        print("To install missing packages, run:")
        print("  pip install -r {}".format(os.path.join(os.path.dirname(__file__), 'requirements.txt')))
        return False
    
    return True

def check_pipeline_files():
    """Check if PairTCR pipeline files exist"""
    web_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(web_dir)
    
    required_files = [
        os.path.join(parent_dir, 'scripts', 'pipeline_runner_base.py'),
        os.path.join(parent_dir, 'scripts', 'mixcr.jar'),
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print("Missing required PairTCR pipeline files:")
        for file_path in missing_files:
            print("  - {}".format(file_path))
        return False
    
    return True

def main():
    parser = argparse.ArgumentParser(description='Start PairTCR Web Server')
    parser.add_argument('--host', default='0.0.0.0', help='Host to bind to (default: 0.0.0.0)')
    parser.add_argument('--port', type=int, default=5000, help='Port to bind to (default: 5000)')
    parser.add_argument('--debug', action='store_true', help='Run in debug mode')
    parser.add_argument('--check-only', action='store_true', help='Only check dependencies, do not start server')
    
    args = parser.parse_args()
    
    print("PairTCR Web Server Startup")
    print("=" * 50)
    
    # Check Python version
    print("Checking Python version...")
    if not check_python_version():
        sys.exit(1)
    print("✓ Python version OK")
    
    # Check dependencies
    print("Checking dependencies...")
    if not check_dependencies():
        sys.exit(1)
    print("✓ All dependencies installed")
    
    # Check pipeline files
    print("Checking PairTCR pipeline files...")
    if not check_pipeline_files():
        sys.exit(1)
    print("✓ Pipeline files found")
    
    if args.check_only:
        print("\n✓ All checks passed!")
        return
    
    # Start the server
    print("\nStarting PairTCR Web Server...")
    print("Host: {}".format(args.host))
    print("Port: {}".format(args.port))
    print("Debug: {}".format(args.debug))
    print("=" * 50)
    
    try:
        # Change to web directory
        web_dir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(web_dir)
        
        # Import and run the Flask app
        from app import app, socketio
        
        print("\n🚀 Server starting...")
        print("📱 Access the application at: http://localhost:{}".format(args.port))
        print("🔗 Or from another machine: http://{}:{}".format(args.host, args.port))
        print("\nPress Ctrl+C to stop the server")
        print("=" * 50)
        
        socketio.run(app, host=args.host, port=args.port, debug=args.debug)
        
    except KeyboardInterrupt:
        print("\n\n🛑 Server stopped by user")
    except Exception as e:
        print("\n❌ Error starting server: {}".format(e))
        sys.exit(1)

if __name__ == '__main__':
    main() 