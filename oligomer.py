#!/usr/bin/python

import re

# you can make a regular expression for the match
# leave out the ^ if you don't want to only match it
# at the start of the line
# (The ` characters turn int into strings so you can concatenate them in filenames and other strings, etc.)  
target_pattern=re.compile('Oligomer2') 

outfile=open('Oligomer2.txt', 'w')

#counter=0

for n in range(1,1000):
  filein = 'run_' + str(n) + '.txt'
  f=open(filein,'r')

  size=-1

  for line in f:
    found=re.search(target_pattern, line)
    if found:
      size=21
      #counter+=1
      # size will now be the flag/counter to capture the current line while looping through
      # and we'll reduce by 1 later on each time a line is sent out
    # this will spit out the current line with the marker, but if you don't want that
    # then put this block of code in front of the 'if found' part
    if size==0:
      # get one last line (delete this part if you move the block to not have the marker
      # line included and you then the counter would be offset by 1)
      outfile.write(line)
      size=size-1
    if size>0:
      outfile.write(line)
      size=size-1
  f.close()

#print counter,
#outfile.write('counter='+`counter`)
outfile.close()
