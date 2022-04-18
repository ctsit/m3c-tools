import datetime
import os.path


class Logger(object):
    def __init__(self, filepath='log.txt', is_verbose=None):
        self.is_verbose = is_verbose

        filename = os.path.basename(filepath)
        date = str(datetime.date.today())
        name = date + '_' + filename
        self.path = os.path.dirname(filepath)
        self.log_path = os.path.join(self.path, name)

        sep = '#' * 10
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        time_message = sep + " " + timestamp + " " + sep
        self.__call__(time_message)
        self.__call__(sep * 4)

    def __call__(self, message):
        with open(self.log_path, "a+") as log:
            log.write(message + "\n")
        if self.is_verbose:
            print(message)

    def get_location(self):
        location = self.path
        return location
