# Get the absolute path of the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if the script directory is already in PATH
if [[ ":$PATH:" == *":$SCRIPT_DIR:"* ]]; then
  echo "The directory is already in PATH. CLIPittyClip is probably already installed."
else
  echo "This installation script will add the current directory to your PATH and give execute permission so CLIPittyClip can be run natively as a command-line program."
  read -n 1 -r -p "Do you want to proceed? (Press any key to continue or 'Ctrl+C' to cancel) "
  echo

  echo "Adding $SCRIPT_DIR to PATH..."
  echo "PATH=$PATH:$SCRIPT_DIR" >> ~/.bashrc
  chmod +x "$SCRIPT_DIR/CLIPittyClip.sh"
  chmod +x "$SCRIPT_DIR/CombinedPeakCalling.sh"
  echo "Done! You can now run CLIPittyClip.sh in directory containing your fastq.gz file."
  echo "Make sure to restart your terminal or run 'source ~/.bashrc'."
fi