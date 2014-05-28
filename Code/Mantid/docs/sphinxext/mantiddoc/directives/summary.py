from base import BaseDirective


class SummaryDirective(BaseDirective):

    """
    Obtains the summary for a given algorithm based on it's name.
    """

    required_arguments, optional_arguments = 0, 0

    def run(self):
        """
        Called by Sphinx when the ..summary:: directive is encountered.
        """
        title = self._make_header("Summary")
        summary = self._get_summary()
        return self._insert_rest(title + summary)

    def _get_summary(self):
        """
        Return the summary for the named algorithm.

        Args:
          algorithm_name (str): The name of the algorithm.
        """
        alg = self._create_mantid_algorithm(self.algorithm_name(), self.algorithm_version())
        return alg.getWikiSummary()


def setup(app):
    """
    Setup the directives when the extension is activated

    Args:
      app: The main Sphinx application object
    """
    app.add_directive('summary', SummaryDirective)
