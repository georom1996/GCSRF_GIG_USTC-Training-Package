#Run the function part in the last box before you test the processes
import os
from gcsrf_plot_subroutines import *
calculate_path="./calculate_test/"
sub_fd_lst=sorted(os.listdir(calculate_path))
GCSRF_plotSyn(calculate_path,sub_fd_lst,"syn","(a)","(b)","(c)","(d)","(e)","syn-SM-d75")