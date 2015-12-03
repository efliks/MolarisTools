#-------------------------------------------------------------------------------
# . File      : Utilities.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------

def TokenizeLine (line, converters=None, separator=None, reverse=False):
    """Tokenize a line with optional converters and separator."""
    tokens = None
    if line is not None:
        if separator is None:
            tokens = line.split ()
        else:
            tokens = line.split (separator)
        # . The last token becomes first
        if reverse:
            tokens.reverse ()
        if converters is not None:
            ntokens     = len (tokens)
            nconverters = len (converters)
            # . Truncate tokens if there is more of them than converters
            if ntokens > nconverters:
                tokens = tokens[:nconverters]
            # . Add "empty" tokens if there is fewer of them than converters
            if ntokens < nconverters:
                tokens.extend ([None] * (nconverters - ntokens))
            # . Do the conversion
            for (i, (token, converter)) in enumerate (zip (tokens, converters)):
                if (token is None) or (converter is None):
                    new = token
                else:
                    try:
                        new = converter (token)
                    except:
                        raise exceptions.StandardError ("Unable to convert token " + repr (i) + ".", True)
                tokens[i] = new
    return tokens

def WriteData (data, filename, append=False):
        output = open (filename, "a" if append else "w")
        output.writelines (data)
        output.close ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
