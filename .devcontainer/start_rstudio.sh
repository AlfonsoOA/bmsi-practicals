#!/bin/bash
# Runs every time the codespace window opens: starts RStudio Server (as the user,
# never with sudo) and prints the address to open it. The address is the
# forwarded one (https://...-8787.app.github.dev), which works in every browser;
# "localhost:8787" only works with VS Code installed on the computer.
rstudio-server start > /dev/null 2>&1 || true
URL="https://${CODESPACE_NAME}-8787.${GITHUB_CODESPACES_PORT_FORWARDING_DOMAIN:-app.github.dev}"
echo "$URL" > ~/RSTUDIO_ADDRESS.txt
printf '\n=====================================================================\n'
printf '  RStudio is ready. Open it here (Ctrl + click on the address):\n\n'
printf '  %s\n' "$URL"
printf '=====================================================================\n\n'
