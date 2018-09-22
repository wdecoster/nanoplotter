from base64 import b64encode
from io import BytesIO
from urllib.parse import quote as urlquote


class Plot(object):
    """A Plot object is defined by a path to the output file and the title of the plot."""

    def __init__(self, path, title):
        self.path = path
        self.title = title
        self.fig = None
        self.html = None

    def encode(self):
        if self.html:
            return self.html
        elif self.fig:
            return self.encode2()
        else:
            return self.encode1()

    def encode1(self):
        """Return the base64 encoding of the figure file and insert in html image tag."""
        data_uri = b64encode(open(self.path, 'rb').read()).decode('utf-8').replace('\n', '')
        return '<img src="data:image/png;base64,{0}">'.format(data_uri)

    def encode2(self):
        """Return the base64 encoding of the fig attribute and insert in html image tag."""
        buf = BytesIO()
        self.fig.savefig(buf, format='png', bbox_inches='tight', dpi=100)
        buf.seek(0)
        string = b64encode(buf.read())
        return '<img src="data:image/png;base64,{0}">'.format(urlquote(string))
