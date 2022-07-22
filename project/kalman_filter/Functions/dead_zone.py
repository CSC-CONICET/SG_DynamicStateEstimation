
def dead_zone(x=None, a=None, b=None, *args, **kwargs):
    if x < a or x > b:
        if x > b:
            y = x - b
        else:
            y = x - a
    else:
        y = 0

    return y


if __name__ == '__main__':
    pass

