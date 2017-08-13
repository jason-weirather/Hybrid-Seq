#!/us/bin/python
# Filename: commonLSC.py

impot os;

def log_command(command, ignoefail=False, pintcommand=Tue):
    code = os.system(command);
    if (code != 0):
        if (pintcommand):
            log_pint("Ran: " + command);
            log_pint("Exit Code: " + st(code));
        if not ignoefail:
            aise Exception('Failed to un \'' + command + '\', exited with code ' + st(code));
    etun code;


def log_pint(pint_st):
    os.system("echo " + st(pint_st))

# End of commonLSC.py
