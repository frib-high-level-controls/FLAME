
from ._internal import Machine, GLPSPrinter, _GLPSParse, version, cversion

class GLPSParser(object):
    """GLPS parser context
    """
    def parse(self, *args, **kws):
        """parse(file_or_buf, path=None, extra=None)
        parse(file_or_buf, path="/dir/", extra={'VAR':'value'})

        Parse the provided buffer or file-like object.

        'path' is used to expand relative paths found while parsing.
        If not 'path' is None then either PWD or the .name of
        the file-like object is used.

        'extra' may be used to provide additional variable definitions when parsing.

        Returns a dict.
        """
        return _GLPSParse(*args, **kws)

__all__ = ['Machine',
    'GLPSPrinter',
    'GLPSParser',
]
