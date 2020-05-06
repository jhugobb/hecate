import os
import subprocess
import sys
import yaml

# ./run.py hecate_bash_path folder_with_models configfile
def main():
  if len(sys.argv) != 3: 
    print("Usage: ./run.py folder_with_models_hecs configfile")
    return -1

  args = sys.argv
  folder = args[1]
  configfile = args[2]

  with open(configfile, 'r') as stream:
    try:
      config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
      print(exc)

  outputfile = config["output"]
  outputfile_stream = open(outputfile, 'w')
  use_7z = config["7z"]
  use_gzip = config["gzip"]

  for filename in os.listdir(folder):
    if filename.endswith(".hec"):
      print("Running " + filename)

        # cmd_str = hecate_bash + " " + os.path.relpath(folder+"/"+filename) + " " + str(grid_size) + " " + str(voxelization_tech) + " " + str(calc_BW*1) + " " + \
        #           str(threshold_ray) + " " + str(writePNG*1) + " " + str(writePLY*1) + " " + str(writeCSV*1) + " " + str(writeHEC*1)
        
      if use_7z:
        cmd_str = "7z a " + folder+"/"+filename[:-4] + ".7z " + os.path.relpath(folder+"/"+filename)
        if subprocess.call(cmd_str,stdout=outputfile_stream, shell=True) != 0:
          print("Error when running " + filename +"!")
          exit(2)
      
      if use_gzip:
        cmd_str = "gzip -ck " + os.path.relpath(folder+"/"+filename) + " > " + folder+"/"+filename[:-4] + ".gz " 
        if subprocess.call(cmd_str,stdout=outputfile_stream, shell=True) != 0:
          print("Error when running " + filename +"!")
          exit(2)
    else:
      continue
  
  outputfile_stream.close()

if __name__ == "__main__":
  main()