import glob
import os
import time
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
      total_time_7z = 0.0
      total_time_gz = 0.0
      for candidate_file in os.listdir(folder+"/"+filename):
        if candidate_file.endswith(".hec"):
          # print("Running " + filename+"/"+candidate_file)

          start_7z = time.clock_gettime(time.CLOCK_REALTIME)
          if use_7z:
            cmd_str = "7z a " + folder+"/"+filename+"/"+candidate_file[:-4] + ".7z " + os.path.relpath(folder+"/"+filename+"/"+candidate_file)
            if subprocess.call(cmd_str,stdout=outputfile_stream, shell=True) != 0:
              print("Error when running " + candidate_file +"!")
              exit(2)
          end_7z = time.clock_gettime(time.CLOCK_REALTIME)
          total_time_7z += end_7z - start_7z

          start_gz = time.clock_gettime(time.CLOCK_REALTIME)
          if use_gzip:
            cmd_str = "gzip -ck " + os.path.relpath(folder+"/"+filename+"/"+candidate_file) + " > " + folder+"/"+filename+"/"+candidate_file[:-4] + ".gz " 
            if subprocess.call(cmd_str,stdout=outputfile_stream, shell=True) != 0:
              print("Error when running " + candidate_file +"!")
              exit(2)
          end_gz = time.clock_gettime(time.CLOCK_REALTIME)
          total_time_gz += end_gz-start_gz
        else:
          continue
      print("7z of " + folder+"/"+filename + " was " + str(total_time_7z))
      print("gz of " + folder+"/"+filename + " was " + str(total_time_gz))

      if (os.listdir(folder+"/"+filename)):
        hec_files = glob.glob(folder+"/"+filename+"/*.hec")
        z7_files = glob.glob(folder+"/"+filename+"/*.7z")
        gz_files = glob.glob(folder+"/"+filename+"/*.gz")
        if subprocess.call("mkdir " + folder+"/"+filename+"/hec " + folder+"/"+filename+"/7z " + folder+"/"+filename+"/gz", stdout=outputfile_stream, shell=True) != 0:
          print("Error while trying to create folders!")
          exit(2)
          "mv " + hec_files + " " + folder+"/"+filename+"/hec"
        if subprocess.call(['mv'] + hec_files + [folder+"/"+filename+"/hec"], stdout=outputfile_stream) != 0:
          print("Error while trying to move hec!")
          exit(2)
        if subprocess.call(['mv'] + z7_files + [folder+"/"+filename+"/7z"], stdout=outputfile_stream) != 0:
          print("Error while trying to move 7z!")
          exit(2)
        if subprocess.call(['mv'] + gz_files + [folder+"/"+filename+"/gz"], stdout=outputfile_stream) != 0:
          print("Error while trying to move gz!")
          exit(2)
  
  outputfile_stream.close()

if __name__ == "__main__":
  main()