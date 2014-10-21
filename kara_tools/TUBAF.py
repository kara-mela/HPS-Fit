def blau(style='TUBAF'):
    if style=='TUBAF':
        return (0.000, 0.392, 0.659)
    else:
        return 'orange'

def gruen(style='TUBAF'):
    if style=='TUBAF':
        return (0.102, 0.588, 0.168)
    else:
        return (170./255, 218./255, 0./255) # bright green
        
def rot(style='TUBAF'):
    if style=='TUBAF':
        return (0.710, 0.071, 0.243)
    else:
        return (0./255, 162./255, 232./255) # bright blue
        
def cyan(style='TUBAF'):
    if style=='TUBAF':
        return (0.137, 0.729, 0.886)
    else:
        return 'red'
        
def orange(style='TUBAF'):
    if style=='TUBAF':
        return (0.878, 0.525, 0.112)
    else:
        return 'blue'

def gelb(style='TUBAF'):
    if style=='TUBAF':
        return (216./255, 188./255, 61./255)
    else:
        return (0.196, 0.804, 0.196) # forest green

def width(style='TUBAF'):
    if style=='TUBAF':
        return 1.5
    else:
        return 1
        
def name(style='TUBAF'):
    if style=='TUBAF':
        return 'TUBAF'
    else:
        return 'standard'