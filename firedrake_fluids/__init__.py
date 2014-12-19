import logging

log = logging.getLogger(__name__)
_handler = logging.StreamHandler()
log.addHandler(_handler)
_handler.setFormatter(logging.Formatter('%(name)s %(levelname)s: %(message)s'))
del(_handler)
log.setLevel(logging.DEBUG)
