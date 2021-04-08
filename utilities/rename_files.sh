#! /bin/bash
find . -depth -name '*sound_velocity*' -execdir bash -c 'mv -i "$1" "${1//sound_velocity/soundVelocity}"' bash {} \;