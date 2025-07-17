import builtins

def verbose_print(f):
    def wrap(*args,**kwargs):
        # Deal with verbose keyword
        if 'verbose' in kwargs.keys():
            verbose = kwargs['verbose']
            kwargs.pop("verbose")
        else:
            verbose = False
        
        # redefine the print function to only run when we trigger verbose
        p = builtins.print
        def log(*print_args,**print_kwargs):
            if verbose:
                p(*print_args,**print_kwargs)
            else:
                pass
        
        # Run function
        builtins.print = log
        res = f(*args,**kwargs)
        
        # Reset print
        builtins.print = p
        
        return res
    return wrap

