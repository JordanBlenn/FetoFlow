file_parsing_utils
===================

Utilities for reading the simple ipnode/ipelem/ipfiel file formats used by
the project. The functions here are intentionally small and strict — they
perform basic validation of file names and convert the 1-based indexing of
the input files to the 0-based indexing used internally.

Functions
---------

read_nodes(filename)
    Parse a ``.ipnode`` file and return a dictionary mapping node ids (0-based)
    to their (x, y, z) coordinates. Performs basic filename validation and
    converts coordinates to floats.
    Parameters
    ----------
    filename : str
        Path to a ``.ipnode`` file. The function validates the extension
        and expects the format used by the project's mesh exporter.

    Returns
    -------
    dict
        Mapping ``node_index -> [x, y, z]`` with 0-based node indices.

    Examples
    --------
    >>> nodes = read_nodes('placenta.ipnode')

read_elements(filename)
    Parse a ``.ipelem`` file and return a list of element tuples
    (node_from, node_to) using 0-based node indices.

    Parameters
    ----------
    filename : str
        Path to a ``.ipelem`` file.

    Returns
    -------
    list of tuple
        Each element is ``(node_from_index, node_to_index)`` with 0-based
        indices.

define_fields_from_files(files)
    Read one or more ``.ipfiel`` files that encode per-element/per-node
    fields (for example a radius or resistance field). The input is a
    dictionary mapping field names to filenames and the return value is a
    dict mapping field names to dictionaries of element indices -> values.

    Parameters
    ----------
    files : dict
        Dictionary mapping field names (``str``) to filenames (``str``)
        of corresponding ``.ipfiel`` files.

    Returns
    -------
    dict
        Mapping field name to a dict of element indices -> numeric value.

    Notes
    -----
    The function performs minimal validation of the ipfiel format and
    will raise :class:`TypeError` for unexpected extensions.

Notes and limitations
---------------------

- The parser expects the original simple text formats used in the
  project. If you need to support other mesh/file formats consider
  converting them upstream or adding another helper.
- The functions validate file extensions and will raise :class:`TypeError`
  for unexpected extensions.

Examples
--------

.. code-block:: python

    nodes = read_nodes('placenta.ipnode')
    elements = read_elements('placenta.ipelem')
    fields = define_fields_from_files({'radius':'radius.ipfiel'})

API reference
-------------

.. automodule:: FetoFlow.file_parsing_utils
    :members:
    :undoc-members:
    :show-inheritance:


read_nodes
    Parameters
    ----------
    filename : str
        Path to a ``.ipnode`` file. The function validates the extension
        and expects the format used by the project's mesh exporter.

    Returns
    -------
    dict
        Mapping ``node_index -> [x, y, z]`` with 0-based node indices.

read_elements
    Parameters
    ----------
    filename : str
        Path to a ``.ipelem`` file. The function validates the extension
        and returns a list of tuples describing element connectivity.

    Returns
    -------
    list of tuple
        Each element is ``(node_from_index, node_to_index)`` with 0-based
        indices.

define_fields_from_files
    Parameters
    ----------
    files : dict
        Dictionary mapping field names (``str``) to filenames (``str``)
        of corresponding ``.ipfiel`` files.

    Returns
    -------
    dict
        Mapping field name to a dict of element indices -> numeric value.

Notes
-----
These functions are intentionally strict; if you need to import other
mesh formats create a lightweight converter that outputs the simple
``.ipnode``/``.ipelem`` format and call the helpers here. See
:mod:`FetoFlow.geometry_utils` for how parsed nodes/elements are used to
construct the graph.

Examples
--------

.. code-block:: python

    nodes = read_nodes('input/placenta.ipnode')
    elements = read_elements('input/placenta.ipelem')
    fields = define_fields_from_files({'radius':'input/radius.ipfiel'})

