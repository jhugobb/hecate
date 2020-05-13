import os
import subprocess
import sys
import yaml

# ./run.py hecate_bash_path folder_with_models configfile
def main():
  if len(sys.argv) != 3: 
    print("Usage: ./run.py folder_with_dirs_hecs configfile")
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
    # print ("I entered " + str(os.path.isdir(os.path.relpath(filename))))
    if os.path.isdir(folder+"/"+filename):
      for candidate_file in os.listdir(folder+"/"+filename):
        if candidate_file.endswith(".hec"):
          print("Running " + filename+"/"+candidate_file)

          if use_7z:
            cmd_str = "7z a " + folder+"/"+filename+"/"+candidate_file[:-4] + ".7z " + os.path.relpath(folder+"/"+filename+"/"+candidate_file)
            if subprocess.call(cmd_str,stdout=outputfile_stream, shell=True) != 0:
              print("Error when running " + candidate_file +"!")
              exit(2)
          
          if use_gzip:
            cmd_str = "gzip -ck " + os.path.relpath(folder+"/"+filename+"/"+candidate_file) + " > " + folder+"/"+filename+"/"+candidate_file[:-4] + ".gz " 
            if subprocess.call(cmd_str,stdout=outputfile_stream, shell=True) != 0:
              print("Error when running " + candidate_file +"!")
              exit(2)
        else:
          continue
  
  outputfile_stream.close()

if __name__ == "__main__":
  main()