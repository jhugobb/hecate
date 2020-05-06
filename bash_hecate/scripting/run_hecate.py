import os
import subprocess
import sys
import yaml

# ./run.py hecate_bash_path folder_with_models configfile
def main():
  if len(sys.argv) != 4: 
    print("Usage: ./run_hecate.py hecate_bash_path folder_with_models configfile")
    return -1

  args = sys.argv
  hecate_bash = args[1]
  folder = args[2]
  configfile = args[3]

  with open(configfile, 'r') as stream:
    try:
      config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
      print(exc)

  outputfile = config["output"]
  outputfile_stream = open(outputfile, 'w')
  sizes = config["grid_size"]
  calc_BW = config["calculate_BW"]
  voxelization_tech = config["voxelization_tech"]
  threshold_ray = config["thresh_raycasting"]
  writePNG = config["png"]
  writePLY = config["ply"]
  writeCSV = config["csv"]
  writeHEC = config["hec"]

  for filename in os.listdir(folder):
    if filename.endswith(".ply"):
      print("Running " + filename)
      for grid_size in sizes:
        print("Running for " + str(grid_size))
        cmd_str = hecate_bash + " " + os.path.relpath(folder+"/"+filename) + " " + str(grid_size) + " " + str(voxelization_tech) + " " + str(calc_BW*1) + " " + \
                  str(threshold_ray) + " " + str(writePNG*1) + " " + str(writePLY*1) + " " + str(writeCSV*1) + " " + str(writeHEC*1)
        if subprocess.call(cmd_str,stdout=outputfile_stream, shell=True) != 0:
          print("Error when running " + filename +"!")
          exit(2)

    else:
      continue
  
  outputfile_stream.close()
if __name__ == "__main__":
  main()