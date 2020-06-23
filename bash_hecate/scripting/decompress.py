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

  for filename in os.listdir(folder):
    # print ("I entered " + str(os.path.isdir(os.path.relpath(filename))))
    if os.path.isdir(folder+"/"+filename):
      for candidate_file in os.listdir(folder+"/"+filename):
        if os.listdir(folder+"/"+filename+"/"+candidate_file) and candidate_file == "7z":
          total_time_7z = 0.0
          print("Running 7z for " + filename)
          for comp_file in os.listdir(folder+"/"+filename+"/"+candidate_file):
            if comp_file.endswith(".7z"):
              start_7z = time.clock_gettime(time.CLOCK_REALTIME)
              cmd_str = "7z e " + folder+"/"+filename+"/"+candidate_file+"/"+comp_file+ " -o" + os.path.relpath(folder+"/"+filename+"/"+candidate_file)
              if subprocess.call(cmd_str,stdout=outputfile_stream, shell=True) != 0:
                print("Error when decompressing " + comp_file +"!")
                exit(2)
              end_7z = time.clock_gettime(time.CLOCK_REALTIME)
              total_time_7z += end_7z - start_7z
          print("7z of " + folder+"/"+filename + " was " + str(total_time_7z))
          hec_files = glob.glob(folder+"/"+filename+"/"+candidate_file+"/*.hec")
          png_files = glob.glob(folder+"/"+filename+"/"+candidate_file+"/*.png")
          if len(hec_files) != 0 and subprocess.call(['rm'] + hec_files, stdout=outputfile_stream) != 0:
            print("Error while trying to remove hec!")
            exit(2)
          if len(png_files) != 0 and subprocess.call(['rm'] + png_files, stdout=outputfile_stream) != 0:
            print("Error while trying to remove png!")
            exit(2)
        elif os.listdir(folder+"/"+filename+"/"+candidate_file) and candidate_file == "gz":
          print("Running gz for " + filename)
          total_time_gz = 0.0
          for comp_file in os.listdir(folder+"/"+filename+"/"+candidate_file):
            if comp_file.endswith(".gz"):
              start_gz = time.clock_gettime(time.CLOCK_REALTIME)
              cmd_str = "gzip -kdN " + folder+"/"+filename+"/"+candidate_file+"/"+comp_file
              if subprocess.call(cmd_str,stdout=outputfile_stream, shell=True) != 0:
                print("Error when decompressing " + comp_file +"!")
                exit(2)
              end_gz = time.clock_gettime(time.CLOCK_REALTIME)
              total_time_gz += end_gz - start_gz
          print("gz of " + folder+"/"+filename + " was " + str(total_time_gz))
          hec_files = glob.glob(folder+"/"+filename+"/"+candidate_file+"/*.hec")
          png_files = glob.glob(folder+"/"+filename+"/"+candidate_file+"/*.png")
          if len(hec_files) != 0 and subprocess.call(['rm'] + hec_files, stdout=outputfile_stream) != 0:
            print("Error while trying to remove hec!")
            exit(2)
          if len(png_files) != 0 and subprocess.call(['rm'] + png_files, stdout=outputfile_stream) != 0:
            print("Error while trying to remove png!")
            exit(2)
        else:
          continue
  outputfile_stream.close()

if __name__ == "__main__":
  main()