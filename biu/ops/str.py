def safe(s):
    """
    Make a string a valid identifier
    """
    sl = 'abcdefghijklmnopqrstuvwxyz' + 'abcdefghijklmnopqrstuvwxyz'.upper() + '_1234567890'
    whitespace = ' \n\t\r'
    s = ''.join([ l if l in sl else '_' if l in whitespace else '' for l in s]).strip()
    if s[0].isnumeric():
        s = 'X' + s
    #fi
    return s
#edef