import sys
import os
import subprocess

def run_one_test(filename, output):
    subprocess.check_call(["python", filename],
                          stderr=subprocess.STDOUT,
                          stdout=output,
                          timeout=10)

# I don't use unittest runner because I want to run each test in a separate
# subprocess; this is because OPAL instance is global
def main():
    passes, fails = [], []
    this_dir = os.path.split(__file__)[0]
    here = os.getcwd()
    if this_dir == "":
        this_dir = "."
    fout = open(os.path.join(this_dir, "test.log"), "w")
    for dirpath, dirnames, filename_list in os.walk(this_dir):
        os.chdir(dirpath)
        for filename in filename_list:
            if filename[:5] == "test_":
                try:
                    run_one_test(filename, fout)
                    passes.append(filename)
                except subprocess.CalledProcessError:
                    print(filename, "fails")
                    fails.append(filename)
        os.chdir(here)
    output_string = "Passes:\n"
    for fname in passes:
        output_string += "    "+fname+"\n"
    if len(passes) == 0:
        output_string += "    <None>\n"
    output_string += "Fails:\n"
    for fname in fails:
        output_string += "    "+fname+"\n"
    if len(fails) == 0:
        output_string += "    <None>\n"
    output_string += str(len(passes))+" passed and "+str(len(fails))+" failed"
    print(output_string)
    print(output_string, file = fout)
    return len(fails) > 0


if __name__ == "__main__":
    fails = main()
    sys.exit(fails)
    