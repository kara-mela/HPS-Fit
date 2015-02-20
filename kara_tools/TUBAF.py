import pylab as pl

def width(style='TUBAF'):
    if style=='TUBAF':
        return 2
    else:
        return 1.5
        
def name(style='TUBAF'):
    if style=='TUBAF':
        return 'TUBAF'
    else:
        return 'standard'

def color(style='TUBAF'):
    """
    style==TUBAF:
    b - blue, g - green, r - red, c - cyan, o - orange, y - yellow
    else:
    b - orange, g - brightgreen, r - bright blue, c - red, o - blue, y - borest green
    """
    if style == 'TUBAF':
        colors = dict({'b' : (0.000, 0.392, 0.659),
                       'g' : (0.102, 0.588, 0.168),
                       'r' : (0.710, 0.071, 0.243),
                       'c' : (0.137, 0.729, 0.886),
                       'o' : (0.878, 0.525, 0.112),
                       'y' : (0.847, 0.737, 0.239)})
    else:
        colors = dict({'b' : 'orange',
                       'g' : (0.55, 0.85, 0.000), # (0.667, 0.855, 0.000), # bright green
                       'r' : (0.000, 0.635, 0.910), # bright blue
                       'c' : 'red',
                       'o' : 'blue',
                       'y' : (0.196, 0.804, 0.196) # forest green
                       })
    return colors
    
def cl(style='TUBAF'):
    my_c = color(style)
    cl = [my_c["r"], my_c["o"], my_c["y"], my_c["g"], my_c["b"], my_c["c"]]
    return cl

def help_lines(pos, limits=[-1, 105]):
    try:
        if type(pos) == type(1):
            pl.plot([pos, pos], [limits[0], limits[1]], color='gray', lw=1.5, linestyle='--')
        else:
            for i in pos:
                pl.plot([i, i], [limits[0], limits[1]], color='gray', lw=1.5, linestyle='--')
    except ValueError:
        print("Invalid data type for %s"%pos)