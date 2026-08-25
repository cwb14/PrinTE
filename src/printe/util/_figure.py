"""Shared figure options.

Journals strip figure titles and carry the text in the caption instead, so PrinTE's
plots are untitled by default. Passing --title on its own puts the script's own
descriptive title back for reading on screen; --title "text" overrides it.
"""


def add_title_args(parser):
    parser.add_argument(
        "--title", nargs="?", const=True, default=None, metavar="TEXT",
        help="Add a figure title. Bare --title uses this plot's default text; "
             "--title TEXT sets your own. Omitted entirely by default, for publication.",
    )
    return parser


def resolve_title(args, default):
    """Return the title to draw, or None.

    default is the descriptive string the script would have hardcoded.
    """
    chosen = getattr(args, "title", None)
    if chosen is None:
        return None
    return default if chosen is True else chosen
