#!/bin/bash
# Runs every time the codespace window opens: starts RStudio Server (as the user,
# never with sudo) and prints the address to open it. The address is the
# forwarded one (https://...-8787.app.github.dev), which works in every browser;
# "localhost:8787" only works with VS Code installed on the computer.
rstudio-server start > /dev/null 2>&1 || true
# Wait until RStudio really answers (up to 2 minutes), so that "ready" means ready
printf 'Starting RStudio'
READY=no
for i in $(seq 1 60); do
  if curl -s -o /dev/null --max-time 2 http://localhost:8787; then READY=yes; break; fi
  printf '.'; sleep 2
done
echo
URL="https://${CODESPACE_NAME}-8787.${GITHUB_CODESPACES_PORT_FORWARDING_DOMAIN:-app.github.dev}"
echo "$URL" > ~/RSTUDIO_ADDRESS.txt
printf '\n=====================================================================\n'
if [ "$READY" = yes ]; then
  printf '  RStudio is ready. Open it here (Ctrl + click on the address):\n\n'
else
  printf '  RStudio is not answering yet. Type  rstudio-server start  here, wait a\n'
  printf '  few seconds and open this address (Ctrl + click):\n\n'
fi
printf '  %s\n' "$URL"
printf '=====================================================================\n\n'
