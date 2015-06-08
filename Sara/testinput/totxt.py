# Overwrite density_project.data with any of the image files.
# From one directory below:
#  python testinput/totext.py testinput/planet.png > density_project.data

import subprocess
import sys

from string import maketrans

#print sys.argv

command = ['convert', sys.argv[1], 'txt:-']
#print ' '.join(command)
conversion = subprocess.Popen(command, stdout=subprocess.PIPE)

(output, error) = conversion.communicate()

lines = output.split('\n')

empty_trans = maketrans("", "")
for line in lines:
  if line.startswith('#') or not line:
    continue

  rgb = line.split("(")[1].split(")")[0]
  colors = rgb.split(',')
  color = (float(colors[0]) / 255.0 + float(colors[1]) / 255.0 +
      float(colors[2]) / 255.0) / 3.0
  color = 0.03 * (1.0 - color)
  print '   %e ' % color
