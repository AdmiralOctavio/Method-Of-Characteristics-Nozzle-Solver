import subprocess
import sys
import os

def start_app():
    # Get the path to the app.py file relative to this script
    app_path = os.path.join(os.path.dirname(__file__), "app.py")
    # This runs the 'streamlit run app.py' command as a subprocess
    subprocess.run([sys.executable, "-m", "streamlit", "run", app_path])