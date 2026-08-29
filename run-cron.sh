#!/bin/sh

# Ensure only one instance of /swi/main.py runs at a time
if ! pidof -x "python" >/dev/null; then
    python /swi/main.py
else
    echo "Another instance of main.py is already running. Skipping."
fi
