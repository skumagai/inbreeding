def getnewlinechar(a):
    if a.w:
        return "\r\n"
    else:
        return os.linesep
