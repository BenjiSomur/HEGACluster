def get_nodes(_data):
    _nodes = []
    for _link in _data:
        _aux = _link.split()
        if _aux[0] not in _nodes:
            _nodes.append(_aux[0])
        if _aux[1] not in _nodes:
            _nodes.append(_aux[1])
    return _nodes


def get_links(_data, _ref):
    _links = []
    for _link in _data:
        _aux = _link.split()
        if _aux[0] == _ref:
            _links.append(_aux[1])
    return _links


def binary_table(_links, _nodes):
    _table = []
    for _ref in _links:
        _line = [0] * len(_nodes)
        for _aux in _ref:
            _line[_nodes.index(_aux)] = 1
        _table.append(_line)
    return _table


def weighted_table(_data, _nodes):
    _table = []
    for _node in _nodes:
        _line = [0] * len(_nodes)
        for _link in _data:
            _aux = _link.split()
            if _aux[0] == _node:
                _line[_nodes.index(_aux[1])] = int(_aux[2])
        _table.append(_line)
    return _table


def create_table(_data):
    _nodes = get_nodes(_data)
    if len(_data[0].split()) == 2:
        _links = [get_links(_data, _node) for _node in _nodes]
        return binary_table(_links, _nodes)
    elif len(_data[0].split()) == 3:
        return weighted_table(_data, _nodes)
    else:
        return None


if __name__ == '__main__':
    print('Parser Module')
