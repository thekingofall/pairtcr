#!/usr/bin/env python3
"""
PairTCR Web Interface Launcher
Quick launcher for the PairTCR web application from the project root
"""

import os
import sys
import subprocess
from pathlib import Path

def main():
    """Launch PairTCR web application"""
    
    # Get the project root directory
    project_root = Path(__file__).parent.absolute()
    web_dir = project_root / "web"
    
    # Check if web directory exists
    if not web_dir.exists():
        print("❌ Error: Web directory not found!")
        print(f"Expected: {web_dir}")
        print("Please ensure the web application is properly installed.")
        sys.exit(1)
    
    # Check if start_server.py exists
    start_script = web_dir / "start_server.py"
    if not start_script.exists():
        print("❌ Error: Web server script not found!")
        print(f"Expected: {start_script}")
        sys.exit(1)
    
    print("🚀 Starting PairTCR Web Interface...")
    print(f"📁 Project root: {project_root}")
    print(f"🌐 Web directory: {web_dir}")
    print("=" * 50)
    
    try:
        # Change to web directory and run the server
        os.chdir(web_dir)
        
        # Forward all command line arguments to the web server
        cmd = [sys.executable, "start_server.py"] + sys.argv[1:]
        
        # Run the web server
        subprocess.run(cmd)
        
    except KeyboardInterrupt:
        print("\n\n🛑 Web server stopped")
    except Exception as e:
        print(f"\n❌ Error starting web server: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main() 