#!/usr/bin/python

# MPB job submitter
# Version 2.0
# JDT 120816
# Edits for rare-earth project by MR

# Usage:
#
# ./(filename).py run [n]
#
# This will create up to n directories and jobs, and submit them to the cluster
# Once cluster scripts complete, this script will be automatically re-run to collect
# the output and combine it in a 

import os, sys, subprocess, pickle, re, time, numpy, math, random

params = {}
job_list = {}


mode = sys.argv[1]


######## utility functions
# this one loops over sets of parameters, calling fn(params) for each combination
# the form of iters is [ [key1, values1], [key2, values2]...]
def doloop(fn,params,iters = [],cond = lambda x: True):

    if len(iters) < 1:
        
        if cond(params) is True:
            fn( params )

    else:
        thislist = iters[0][1]
        thiskey = iters[0][0]

        for val in thislist:
            params[thiskey] = val
            doloop(fn,params,iters[1:],cond)

def make_mpb_header(params):

	output_string = ""
	
	for k in params:
		if script_keyword_list.count(k) == 0:
			output_string += "(define %s %f)\n" % (k,params[k])

	return output_string


# root directory for creating temporary directories to store output
dirroot = os.getcwd()
# set to mob for normal operation; mbpi to assume inversion symmetry
mpbcmd =  r'/home/jt11/sw/bin/mpb'
# queue to use for bsub
queue = "normal_serial"
# number of k-points per job
kpts_per_job = 1

# files to keep
keep_eps = False
keep_fields = False

# params["wy"] = 400
# params["wz"] = 200
# params["hx"] = 150
# params["hy"] = 200
# params["aper"] = 300
params["kx_list"] = numpy.linspace(0.0,0.5,50)
params["runfn"] = """(run-zeven)
(run-zodd)
"""
params['sym'] = ''

# saved parameters will be [kx, ky, kz, save_list[0], save_list[1], etc.]
save_list = ["wy","wz","aspect","aper","ff"]
sym_strings = {"yodd":"YO","zevenyoddfreqs": "ZEYO", "zoddyevenfreqs": "ZOYE","zevenyevenfreqs":"ZEYE","zoddyodd":"ZOYO", "zevenfreqs":"ZE", "freqs":"ALL"}

outfile = "output.dat"

seq = 0
iterlist = [
	["index",[3.55]],
        ["subindex",[2.4]],
	["aspect",[1.0]],
        ["res",[32]],
	["wy",[400]],
	["wz",[220]],
#	["dy",[1000]],
#	["sy",numpy.linspace(-200,200,17)],
	["runfn",["run-yodd"]],  #  ["(run-yodd-zeven)"]], #,"(run-yeven-zodd)"]],
	["kx_num",[8]],
	["kx_idx",[1]],
	["ff",[0.3144]],
        ["aper",numpy.linspace(130,190,61)]
]

def itercond(params):
	#params['aper'] = params['targetwl']/2.0/params['neff']
	#params['targetfreq'] = params['aper'] / params['targetwl']

	kx_list = numpy.linspace(0,0.5,32)

	this_kx = kx_list # numpy.array_split(kx_list,params['kx_num'])[params['kx_idx']]

	this_kx_str = ""

	for kx in this_kx:
		this_kx_str += "(vector3 %f 0 0)" % (kx)

	params['kx_str'] = this_kx_str

	#params["kx_start"] = params["kx_list"][::10][params["kx_idx"]]
	#params["kx_end"] = params["kx_list"][9::10][params["kx_idx"]]
	#params["hy"] = params["wy"]*params['hy_frac']
	return True

def runmpb(params):
	global seq

	seq = seq + 1

	if len(sys.argv) > 2:
		if seq > int(sys.argv[2]):
			return

#	mpb_header = make_mpb_header(params)

#	print mpb_header

	mpb_script ="""
(set! really-output-epsilon? false)

(define wz %(wz)f)
(define wy %(wy)f)
(define aper %(aper)f)
(define index %(index)f)
(define subindex %(subindex)f)


(define hx 99)
(define hy 224)
(define ff %(ff)f)

;(define (rootfun ff)

;(set! hx (* (sqrt ff) aper (/ 1 (sqrt %(aspect)f ) ) ) )
;(set! hy (* (sqrt ff) wy (sqrt %(aspect)f ) ) )


(set! geometry-lattice (make lattice (size 1 16 16)))
(set! resolution %(res)d)

(set! geometry (list 

			(make block
				(center 0 0 (/ (/ wz aper) 2)) (size 1 (/ wy aper) (/ wz aper))
				(material (make dielectric (epsilon (* index index)))))

			(make ellipsoid
				(center 0 0 (/ (/ wz aper) 2)) (size (/ hx aper) (/ hy aper) infinity)
				(material (make dielectric (epsilon 1))))

            (make block                                                                                                              
                (center 0 0 -4) (size 1 16 8)                                                                                     
                (material (make dielectric (epsilon (* subindex subindex)))))

))

(set! num-bands 3)

;(set! k-points (list %(kx_str)s))
(set! k-points (list (vector3 0.5 0 0)))

(%(runfn)s display-zparities)

;(output-at-kpoint (vector3 0.5 0 0) (fix-efield-phase 1) (output-efield 1))
;(output-at-kpoint (vector3 0.5 0 0) (fix-efield-phase 1) (output-efield 2))
;(output-at-kpoint (vector3 0.5 0 0) (fix-efield-phase 1) (output-efield 3))


;;;;; These commands are to 
;;;;; (find-k p omega band-min band-max kdir tol kmag-guess kmag-min kmag-max [band-func...])
;(print "zevenkvals\n")
;(print "zoddkvals\n")
;(print "kvals\n")
;(print "zparity\n")
;(find-k %(sym)s 1.0 1 4 (vector3 1 0 0) 1e-4 0.4 0.01 3.46 display-zparities) ; fix-efield-phase output-efield output-hfield)

""" % params

	#print mpb_script

	# create a temporary directory to hold the script and output
	params['dir_str'] = "mpbpy.%d" % random.randint(1,1e6)
	this_dir = dirroot+ "//" + params['dir_str']
	os.mkdir(this_dir)
	os.chdir(this_dir)

	slurm_script = """#!/bin/bash                                                                                                                            
                                                                                                                                                                 
#SBATCH -n 1 #Number of cores                                                                                                                                    
#SBATCH -t 1440 #Runtime in minutes                                                                                                                                

#SBATCH --mem=5000 #Memory per node in MB (see also --mem-per-cpu)                                                                                              
                                                                                                                                                                 
%(mpbcmd)s mpb_script.mpb > output.txt 2> errors.txt                                                                                                             
python "../%(myname)s" fetch "%(this_dir)s"                                                                                                                      
                                                                                                                                                                 
""" % {'mpbcmd': mpbcmd, 'myname': sys.argv[0], 'this_dir': this_dir}

	# write the script file
	f = open(this_dir + "//mpb_script.mpb","w")
	#f.write(mpb_header)
	f.write(mpb_script)
	f.close()

	# write a script file that will run MPB and then run this script again to get the results                                                                
	exec_script = this_dir + "//mpb_script.sh"
	f = open(exec_script,"w")
	f.write(slurm_script)
	#f.write("%s mpb_script.mpb > mpboutput\n" %mpbcmd)                                                                                                      
	#f.write("python \"../%s\" fetch \"%s\"\n" % (sys.argv[0],this_dir))                                                                                     
	#f.write("h5tovtk *.h5")                                                                                                                                     
	f.close()
	os.chmod(exec_script, 0755)

	#mpb_cmd_str = "echo hello"                                                                                                                              
	#cmd_str = "rbsub -R \"rusage[mem=10000]\" -q %s -o output -e error ./mpb_script.sh" % (queue)                                                           
	cmd_str = "sbatch mpb_script.sh"
	# # now call the LSF queue function to submit the job                                                                                                    
	proc = subprocess.Popen(cmd_str,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)                                                              
	proc.wait()    

	# # record the job ID number for later
	# outstring = proc.stdout.read()
	# #print "Got string: %s" % outstring
	# #output looks like this:
	# # Job <17423602> is submitted to queue <short_serial>.
	# outstringre = re.search("<([0-9]+)>",outstring)
	# jobnum = int(outstringre.group()[1:-1])

	# #print "Got num: %d" % jobnum
	
	# # save job in list; jobnum is assumed to be globally unique so we don't have to worry about conflicts
	# job_list[jobnum] = this_dir
	
	# also write the params to be retrieved later
	f = open(this_dir + "//mpb_script_params.mpb","w")
	pickle.dump(params,f)
	f.close()

	os.chdir(dirroot)

	
# def wait_for_jobs():
# """ take the list of job numbers and poll using bjobs to see if they are finished. This function returns when all the jobs
# have finished"""
# 
# 	alldone=False
# 	while alldone is False:
# 		alldone = True
# 		n_running = 0
# 		n_done = 0
# 		time.sleep(5.0)
# 
# 		for j in job_list:
# 			proc = subprocess.Popen(["bjobs","%d" % j],stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False)
# 			proc.wait()
# 			outstring = proc.stdout.read()
# 			#print "Got string: %s" % outstring
# 
# 			outres = re.match(".*DONE.*",outstring,re.DOTALL)
# 
# 			if len(outstring) > 50 and (outres is None):
# 				n_running = n_running + 1
# 				alldone = False
# 			else:
# 				n_done = n_done + 1
# 
# 		print "%d of %d (%d) processes finished..." % (n_done, n_running+n_done, len(job_list))
#

def retrieve_job():
	""" now that all the jobs are finished, go into each directory and harvest the output, and clean up whatever
	files are necessary. Expect dir to be the directory we are working in """

	os.chdir(dir)

	f = open("output.txt",'r')
	data = f.read()
	f.close()
	
	# this will populate params
	#print params
	f = open("mpb_script_params.mpb",'r')
	params = pickle.load(f)
	f.close()
	
	#print params
	
	f = open(r'..//' + outfile,'a')

	# find all lines containing frequency information
	for iS,iSv in sym_strings.iteritems():
		res = re.findall("^%s.*" % (iS),data,re.MULTILINE)[1:]
		res_split = [s.split(',')[1:] for s in res]
		res_split_num = [ [float(ss) for ss in s] for s in res_split]
		
		print res
		print res_split
		print iS
		print iSv
		
		num_points = len(res_split_num)
		if num_points == 0:
			continue
		num_data = len(res_split_num[0])
		
		# now the first index is the k point, and the second index is
		# 0 k index
		# 1 k1
		# 2 k2
		# 3 k3
		# 4 kmag
		# 5 band1 freq
		# 6 band2 freq, etc
		
		for p in range(num_points):
			f.write("%s,"%iSv)
			
			# write the directory that this data came from, so we can look up files
			f.write("%s,"%params['dir_str'])
			
			for sl in save_list:
				if type(params[sl]) == type('string'):
					f.write("%s,"%params[sl])
				else:
					f.write("%f,"%params[sl])
			
			for pt in range(len(res_split_num[p])):
				f.write("%f,"%res_split_num[p][pt])
			
			f.write("\n")
	
	f.close()


print params
if mode=="run":
	doloop(runmpb,params,iterlist,itercond)
	f = open(dirroot + "//tempdir_list","w")
	pickle.dump(job_list,f)
	f.close()
elif mode=="fetch":
	#dir = dirroot + "//" + sys.argv[2]
	dir = sys.argv[2] # should be complete path
	retrieve_job()


