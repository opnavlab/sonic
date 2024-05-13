import re

def stripFirstLine(raw, name):
    proc = raw
    # Grab text between a `=` and `(` (i.e., a method name)
    regex = re.compile(r'\=(.*?)\(')
    if len(proc) > 0:
        # Execute the regex:
        meth_srch = regex.search(proc[0])
        if meth_srch is not None:
            # Grab the method signature from the regex match:
            meth_sig = meth_srch.group(0)
            # Strip out the method name:
            meth_sig = meth_sig.replace('=', '').replace('(', '').strip()
            # The method name is the final section of the passed-in name:
            meth_name = name.split('.')[-1].strip()
            # If they match, remove the line:
            if meth_sig == meth_name:
                #print(proc[0], ' | ', meth_sig, ' | ', meth_name, ' | ', name)
                del proc[0]

    return proc