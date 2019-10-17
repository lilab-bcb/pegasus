import functools
import time
import gc
import logging

class TimeLogger(object):
    """Log time spent in a function. Can also log an additional custom mesage for debug of info purposes"""
    def __init__(self, custom_message=None, custom_message_type="info"):
        super(TimeLogger, self).__init__()
        self.custom_message = custom_message
        self.custom_message_type = custom_message_type
        self.logger = logging.getLogger("pegasus")
    
    def __call__(self, function):
        @functools.wraps(function)
        def _do(*args, **kwargs):
            start = time.perf_counter()
            res = function(*args, **kwargs)
            end = time.perf_counter()
            
            self.logger.info(
                "Time spent on {0} = {:.2f}s.".format(
                    function.__code__.co_name,
                    end - start
                )
            )

            if self.custom_message is not None:
                log_fct = getattr(self.logger, self.custom_message_type)
                log_fct.info(custom_message)
            return res

        return _do

class GCCollect(object):
    """Force the collection garbage according to binomial distribution"""
    def __init__(self, collection_proba=1):
        super(GCCollect, self).__init__()
        self.collection_proba = collection_proba
    
    def __call__(self, function):
        @functools.wraps(function)
        def _do(*args, **kwargs):
            res = function(*args, **kwargs)
            sample = np.random.binomial(1, self.collection_proba, 1)[0]
            if sample > 0:
                gc.collect()

            return res

        return _do
