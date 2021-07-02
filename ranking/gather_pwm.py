# -*- coding: utf-8 -*-
"""
Created on Oct. 11, 2018

@author: Chun-Ping Yu

class of Positional weight matrix

"""

####
# Usage: python gather_pwm.py collect.meme collect.summary meme_folder
####

import os
class PWM:
	def __init__(self, _id = '', _alength=4, _w=0, _nsites=1, _E = 1.0, _pwm = []):
		self.id = _id
		self.alength = _alength
		self.w = _w
		self.nsites = _nsites
		self.E = _E
		self.pwm = _pwm

	
	def add_by_str(self, _str, _id=None):
		if _id:
			self.id = _str
		else:
			self.id = _id
		self.w = len(_str)
		self.nsites = 1
		self.E = 1.0
		self.pwm = []
		for s in self.id:
			if s == 'A':
				self.pwm.append((1, 0, 0, 0))
			elif s == 'C':
				self.pwm.append((0, 1, 0, 0))
			elif s == 'G':
				self.pwm.append((0, 0, 1, 0))
			elif s == 'T':
				self.pwm.append((0, 0, 0, 1))

	def add_by_value(self, _pwm, _id=None, _nsites=1, _E=1.0):
		if _id:
			self.id = 'Unknown'
		else:
			self.id = _id
		self.w = len(_pwm)
		self.nsites = 1
		self.E = 1.0
		self.pwm = _pwm

	def write(self, _stream):
		_stream.write('MOTIF %s\n\n' % self.id)
		_stream.write('letter-probability matrix: alength= %d w= %d nsites= %d E= %s\n' %
				(self.alength, self.w, self.nsites, self.E))
		#print(self.pwm)
		for row in self.pwm:
			_stream.write('  {:.6f} {:.6f} {:.6f} {:.6f}\n'.format(*row))
		_stream.write('\n')


class SimpleMeme:
	def __init__(self, filename=None):
		if filename:
			self.read(filename)
		else:
			self.background_freq = [0.25, 0.25, 0.25, 0.25]
			self.pwm = []

	def write(self, filename, keep=None):
		with open(filename, 'w') as ofile:
			ofile.write('MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\n')
			value = []
			for c, f in zip(['A', 'C', 'G', 'T'], self.background_freq):
				value.append('{} {:.3f}'.format(c, f))
			ofile.write('%s\n\n' % ' '.join(value))
			if keep:
				for _pwm in self.pwm:
					if _pwm.id in keep:
						_pwm.write(ofile)
			else:
				for _pwm in self.pwm:
					_pwm.write(ofile)

	def read(self, filename):
		self.pwm = []
		for _pwm in self.readPWM(filename):
			self.pwm.append(_pwm)

	def readPWM(self, filename):
		with open(filename) as infile:
			for line in infile:
				if line.startswith('Background'):
					tokens = infile.readline().strip().split()
					self.background_freq = [float(tokens[1]), float(tokens[3]), float(tokens[5]), float(tokens[7])]
					break
			for line in infile:
				if line.startswith('MOTIF'):
					motif = PWM()
					motif.id = line.rstrip().replace('MOTIF ', '')
				elif line.startswith('letter-probability matrix'):
					tokens = line.split()
					motif.alength = int(tokens[3])
					motif.w = int(tokens[5])
					motif.nsites = int(tokens[7])
					motif.E = tokens[9]
					count = 0
					_matrix = []
				elif line.rstrip() == '':
					continue
				else:
					_matrix.append([float(u) for u in line.strip().split()])
					count += 1
					if count == motif.w:
						motif.pwm = _matrix
						yield motif

def rename(name_sites, prefix_name):
	name_sites = sorted(name_sites, key=lambda u: u[1], reverse=True)
	new_name = {}
	for idx, name in enumerate(name_sites):
		new_name[name[0]] = '%s_N%d' % (prefix_name, idx+1)
	return new_name

if __name__ == '__main__':
	import sys, subprocess
	total_pwm = SimpleMeme()
	tmp_file = 'tmp_meme.txt'
	if len(sys.argv) < 4:
		print('Usage: python gather_pwm.py total.meme summary.txt seek_meme_folder')
		exit(1)
	seek_dir = sys.argv[3]
	print('looking for the folder of %s' % seek_dir)
	summary_f = open(sys.argv[2], 'w')
	for root, dirs, files in os.walk('.'):
		if not root.endswith(seek_dir):
			continue
		tmpfile = open(tmp_file, 'w')
		a_file = os.path.join(root, 'meme_out', 'meme.txt')
		if not os.path.exists(a_file):
			continue
		command = ['meme2meme', a_file]
		subprocess.run(command, stdout=tmpfile)
		tmpfile.close()
		pwms = SimpleMeme(tmp_file)
		name_sites = [(a.id, a.nsites) for a in pwms.pwm]
		new_name = rename(name_sites, a_file)
		for a in pwms.pwm:
			a.id = new_name[a.id]
			summary_f.write('%s\t%s\t%d\n' % (a_file, a.id, a.nsites))
		total_pwm.pwm.extend(pwms.pwm)
		print(a_file)
	total_pwm.write(sys.argv[1])
	summary_f.close()
