import os
import signal
import time
from subprocess import Popen, PIPE, STDOUT


def decompress(files_prefix):
    # decompress("tests//exact//exact_")
    #decompress("tests//heuristic//heuristic_")
    for i1 in range(2, 201, 2):
        zeroes = 3 - len(str(i1))
        command = "xz -v -d " + files_prefix
        for i2 in range(zeroes):
            command += '0'
        command += str(i1)
        command += ".gr.xz"
        print(command)
        os.system(command)


def call_cpp(test_prefix, test_id):
    start_time = time.time()
    zeroes = 3 - len(str(test_id))
    input_file = test_prefix
    for j in range(zeroes):
        input_file += '0'
    input_file += str(test_id)
    output_file = input_file + '.tww'
    input_file += '.gr'
    curr_input_test = open(input_file, 'r')
    input_lines = curr_input_test.readlines()
    for line in input_lines:
        tokens = line.strip().split(' ')
        for token in tokens:
            process.stdin.write("{} ".format(token))
        process.stdin.write('\n')
    process.stdin.flush()
    print('\nTest {} is running..\n'.format(input_file), end="")
    output_lines = process.stdout.readlines()
    curr_output_test = open(output_file, 'w')
    for line in output_lines:
        tokens = line.strip().split(' ')
        curr_output_test.writelines('{} {}\n'.format(tokens[0], tokens[1]))
    curr_output_test.close()
    elapsed_time = time.time() - start_time
    print('Test {} has been solved in {} seconds!\n'.format(input_file, elapsed_time), end="")
    print('Solution has {} lines!\n'.format(len(output_lines)), end="")
    print('Solution is verified by the committee checker..')
    os.system("python3 verifier.py {} {}".format(input_file, output_file))
    print("Solution has been verified!\n")


os.system("g++ exact.cpp -o exact ")
#os.system("g++ heuristic.cpp -o heuristic ")
print("Exact solver has been compiled!\n", end="")
print("Heuristic solver has been compiled!\n", end="")
"""
for i in range(1, 11):
    process = Popen('./exact', stdin=PIPE, stdout=PIPE, universal_newlines=True, shell=True, preexec_fn=os.setsid)
    call_cpp(R'tests/tiny/tiny', i)
    os.killpg(os.getpgid(process.pid), signal.SIGTERM)
"""

for i in range(38, 39, 2):
    process = Popen('./exact', stdin=PIPE, stdout=PIPE, universal_newlines=True, shell=True, preexec_fn=os.setsid)
    call_cpp(R'tests/exact/exact_', i)
    os.killpg(os.getpgid(process.pid), signal.SIGTERM)

"""
for i in range(2, 201, 2):
    process = Popen('./heuristic', stdin=PIPE, stdout=PIPE, universal_newlines=True, shell=True, preexec_fn=os.setsid)
    call_cpp(R'tests/heuristic/heuristic_', 200 - i)
    os.killpg(os.getpgid(process.pid), signal.SIGTERM)
"""