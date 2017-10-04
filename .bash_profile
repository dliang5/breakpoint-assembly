# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

#PATH=$PATH:$HOME/bin
PATH=$PATH:/afs/cats.ucsc.edu/users/f/dliang5/bin/sratoolkit/bin:/afs/cats.ucsc.edu/users/f/dliang5/bin/bwa-0.7.15/:/afs/cats.ucsc.edu/users/f/dliang5/bin/reference:/afs/cats.ucsc.edu/users/f/dliang5/bin/samtools-1.3.1/:/afs/cats.ucsc.edu/users/f/dliang5/bin/phrap/:/afs/cats.ucsc.edu/users/f/dliang5/bin/pindel/:~/bin/TE_sites.txt:/afs/cats.ucsc.edu/users/f/dliang5/bin/python2.7/:/afs/cats.ucsc.edu/users/f/dliang5/bin/python3.3
PATH=$PATH+":$HOME/.local/bin:$PATH"
PATH=$PATH+":$HOME/.local1/bin:$PATH" 
#alias python2.7 = "/afs/cats.ucsc.edu/users/f/dliang5/bin/python2.7" 
#alias pytohn3 = "/afs/cats.ucsc.edu/users/f/dliang5/bin/python3.3" 
#export PATH="$HOME/.local/bin:$PATH"
export PATH 
