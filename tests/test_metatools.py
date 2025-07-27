# from unittest.mock import MagicMock, patch

# import pytest

# from bioat.metatools import MetaTools


# @pytest.fixture
# def meta_tools():
#     return MetaTools()


# @patch("bioat.lib.libjgi.JGIOperator")
# @patch("bioat.logger.get_logger")
# def test_JGI_query_happy_path(mock_get_logger, mock_JGIOperator, meta_tools):
#     mock_operator = MagicMock()
#     mock_JGIOperator.return_value = mock_operator

#     meta_tools.JGI_query(query_info="Nemve1", timeout=30, nretry=3)

#     mock_JGIOperator.assert_called_once_with(
#         query_info="Nemve1",
#         xml=None,
#         log_fails=None,
#         nretry=3,
#         timeout=30,
#         regex=None,
#         all_get=False,
#         overwrite_conf=False,
#         filter_files=False,
#         proxy_pool=None,
#         just_query_xml=False,
#         syntax_help=False,
#         usage=False,
#         log_level="DEBUG",
#     )
#     mock_operator.query.assert_called_once()
#     mock_operator.parse_xml.assert_called_once()
#     mock_operator.download.assert_called_once()


# @patch("bioat.lib.libjgi.JGIOperator")
# @patch("bioat.logger.get_logger")
# def test_JGI_query_with_all_get(mock_get_logger, mock_JGIOperator, meta_tools):
#     mock_operator = MagicMock()
#     mock_JGIOperator.return_value = mock_operator

#     meta_tools.JGI_query(query_info="Nemve1", all_get=True)

#     assert mock_operator.all_get is True


# @patch("bioat.lib.libjgi.JGIOperator")
# @patch("bioat.logger.get_logger")
# def test_JGI_query_with_invalid_query_info(
#     mock_get_logger, mock_JGIOperator, meta_tools
# ):
#     with pytest.raises(TypeError):
#         meta_tools.JGI_query(query_info=123)  # Should raise error with non-str input


# @patch("bioat.lib.libjgi.JGIOperator")
# @patch("bioat.logger.get_logger")
# def test_JGI_query_with_no_query_info(mock_get_logger, mock_JGIOperator, meta_tools):
#     mock_operator = MagicMock()
#     mock_JGIOperator.return_value = mock_operator

#     meta_tools.JGI_query(query_info=None)

#     mock_JGIOperator.assert_called_once_with(
#         query_info=None,
#         xml=None,
#         log_fails=None,
#         nretry=4,
#         timeout=60,
#         regex=None,
#         all_get=False,
#         overwrite_conf=False,
#         filter_files=False,
#         proxy_pool=None,
#         just_query_xml=False,
#         syntax_help=False,
#         usage=False,
#         log_level="INFO",
#     )

# @patch("bioat.lib.libjgi.JGIOperator")
# @patch("bioat.logger.get_logger")
# def test_JGI_query_with_negative_timeout(mock_get_logger, mock_JGIOperator, meta_tools):
#     mock_operator = MagicMock()
#     mock_JGIOperator.return_value = mock_operator

#     meta_tools.JGI_query(query_info="Nemve1", timeout=-1)

#     assert mock_operator.timeout == -1  # Should allow negative timeout
