# Get the absolute path of the directory containing this script
SCRIPT_DIR="$(pwd)"

# Check if the script directory is already in PATH
if [[ ":$PATH:" == *":$SCRIPT_DIR:"* ]]; then
  echo "The directory is already in PATH. CLIPittyClip is probably already installed."
else
  echo "This installation script will add the current directory to your PATH and give execute permission so CLIPittyClip programs can be run natively as a command-line program."
  echo "Directory to add: $SCRIPT_DIR"
  echo "Make sure the directory above is the directory containing CLIPittyClip programs. If not, please traverse to the correct directory."
  read -n 1 -r -p "Do you want to proceed? (Press any key to continue or 'Ctrl+C' to cancel) "
  echo

  echo "Adding $SCRIPT_DIR to PATH..."
  echo 'export PATH="$PATH:'$(pwd)'"' >> ~/.bashrc
  chmod +x "$SCRIPT_DIR/CLIPittyClip.sh"
  chmod +x "$SCRIPT_DIR/MAPittyMap.sh"
  chmod +x "$SCRIPT_DIR/PEAKittyPeak.sh"
  echo "Done!"
  echo "Make sure to restart your terminal or run 'source ~/.bashrc'."
fi