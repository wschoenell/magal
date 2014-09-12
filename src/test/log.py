'''
Created on Jul 17, 2012

@author: william

Test the logger
'''

import logging
from magal.core.log import logger

log = logger(__name__)

log.warn('Warn test')
log.info('Info test')
log.error('Err test')
log.debug('dgb test 1')
log.setLevel('DEBUG')
log.debug('dgb test 2')