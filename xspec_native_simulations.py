from subprocess import Popen, PIPE
import time

# Prepare the repeated commands
commands = f"""
parallel leven 2
parallel error 2
mo tbabs*(po+ezdiskbb)
1.0,0
1.7
3.49
0.7
5957.47
fakeit gx339-4_g_low_bgd.pi & gx339-4_g_low.rsp & & y & & test.fak & 1000, 1, 1000
ftgrouppha test.fak backfile=test_bkg.fak test_binned.fak snmin 3
data test_binned.fak
ignore bad
ignore **-2.0 20.0-**
fit
err 1-5
data none
""" * 300 # Repeat 300 times

# Command to start xspec
command = 'xspec'
process = Popen(command, stdin=PIPE, shell=True)

start_time = time.perf_counter()

# Communicate the commands to the process
process.communicate(input=commands.encode())

# Wait for the process to complete if necessary
process.wait()

end_time = time.perf_counter()
total_time = end_time - start_time
print(f"The loop took {total_time} seconds to complete.")