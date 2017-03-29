import sys 

print 'Argument List:', str(sys.argv)
name = str(sys.argv).split(" ") 
write_name = name[1].lstrip("'").rstrip("']") # this is getting the name of the file being written into with actual breakpoints that pass the transponseable element test 
break_point_name = "good_"+name[1].lstrip("'").rstrip("']")
print "this is the writename to the file -?> " + write_name 
print "this is the name of the breakpoint file to check -?> " + break_point_name