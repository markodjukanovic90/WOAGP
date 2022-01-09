import sys, getopt
import random as r

def main(argv):
   
   outputfile = '' # where to save
   m = 0  # params of the instance
   n = 0
   
   try:
      opts, args = getopt.getopt(argv,"ho:m:n:",[ "ofile=", "m=", "n=" ])
   except getopt.GetoptError:
      print('test.py -o <outputfile> -m <m> -n <n>')
      sys.exit(2)
      
   print(opts)
   for opt, arg in opts:
      
      if opt == '-h':
         print('test.py -o <outputfile>')
         sys.exit()
      elif opt in ("-m", "--m"):
         m = arg
      elif opt in ("-n", "--n"):
         n = arg
      elif opt in ("-o", "--ofile", "--m", "--n"):
         outputfile = arg
         
   print('Output file is ', outputfile)
   print('m, n ', m, " ", n)
   
   # generate an instance randomly 
   sets = []
   for i in range(int(m)): # generate i-th set  
       si = set(())
       sets.append(si)
          
   for i in range(int(n)):
       el = r.randint(0, int(m)-1)
       sets[el].add( i + 1 )
       
   # when at least one feasible solution has ensured, add remaining elms randomly in the sets (random lengths of sets are set and filled afterwards)
   for i in range(int(m)):
       up_len = int(n) // 4
       si_len = r.randint(1, up_len) # do not allow possible empty set
       for j in range(si_len):
           el = r.randint(1, int(n) )
           sets[i].add(el)
   
   #for i in range(int(m)):
       #print("Set " + str(i) + ": ")
       #for j in sets[i]:
          #print(str(j), sep=" ")     
   # read into a file
   outfile = "scp_"+ str(n) + "_" + str(m) + ".txt"
   f = open(outfile, "w")
   f.write(n + " " + m)
   f.write("\n")
   for i in range(int(m)):
       si = ""
       for j in sets[i]:
          si +=  str(j) +  " "    
       f.write(si[0:len(si)-1])
       if i < int(m)-1:
          f.write("\n")
       
if __name__ == "__main__":
   main(sys.argv[1:])