import subprocess
import sys 
def main(argv):
    
    ns = [ 10, 50, 100, 500, 1000 ]
    ms = [ [n, n*2, n*5, n*10 ]  for n in ns    ]
    outfile = "test.txt" 
    for n in ns:
         subprocess.call(" python generator_instanci.py -o " + outfile + " -n " + str(n) + " -m " + str(n*30), shell=True)  
       
if __name__ == "__main__":
   main(sys.argv[1:])