#!/bin/bash
PY_SCRIPT_PATH=~/Adafruit_Python_DHT/Adafruit_Python_DHT/examples/
OUTPUT_PATH=~/Adafruit_Python_DHT/output.txt
while true
    do
    OUTPUT=$(sudo ~/Adafruit_Python_DHT/Adafruit_Python_DHT/examples/AdafruitDHT.py 22 17)
#    echo $OUTPUT $(date)
    echo $OUTPUT $(date) >> output.txt
    sleep 60
    done
