import logging
import logging.handlers
import os

import coloredlogs

from .exceptions import BioatInvalidParameterError

__all__ = ["LoggerManager"]


class LoggerManager:
    LOG_FORMAT = "%(asctime)s.%(msecs)03d - [%(name)s] - %(filename)s[line:%(lineno)4d] - %(levelname)+8s: %(message)s"

    def __init__(
        self,
        log_level: str = "ERROR",
        mod_name: str = "bioat.logger",
        cls_name: str | None = None,
        func_name: str | None = None,
    ):
        """
        Initialize the LoggerManager with default log level and module name.

        Args:
            log_level (str, optional): Default logger level. Defaults to "ERROR".
            mod_name (str, optional): Module name. Defaults to "fbt".
            cls_name (str, optional): The name of the class for which the logger is created.
                Defaults to None.
            func_name (str, optional): The name of the function for which the logger is created.
                Defaults to None.
        """
        self.log_level = self._get_log_level(log_level)
        self.mod_name = mod_name
        self.cls_name = cls_name
        self.func_name = func_name
        self.logger = self._create_logger()
        self.logger.propagate = False

    def _get_log_level(self, log_level: str) -> int:
        """Parse and return the logging level."""
        if log_level.upper() not in logging._nameToLevel:
            raise BioatInvalidParameterError(
                f"Invalid log level: {log_level}. "
                "Choose from CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET."
            )
        return logging._nameToLevel[log_level.upper()]

    def _get_logger_name(self) -> str:
        """Construct the logger name based on module, class, and function."""
        if self.cls_name and self.func_name:
            return f"{self.mod_name}.{self.cls_name}.{self.func_name}"
        elif self.cls_name:
            return f"{self.mod_name}.{self.cls_name}"
        elif self.func_name:
            return f"{self.mod_name}.{self.func_name}"
        else:
            return self.mod_name

    def _create_logger(self) -> logging.Logger:
        """Create and configure the logger."""
        name = self._get_logger_name()
        logger = logging.getLogger(name)

        # Console handler with coloredlogs
        coloredlogs.DEFAULT_FIELD_STYLES = dict(
            asctime=dict(color="yellow", bold=True),
            name=dict(color="blue", bold=True),
            filename=dict(color="cyan"),
            lineno=dict(color="green", bold=True),
            levelname=dict(color="magenta", bold=True),
        )
        coloredlogs.DEFAULT_LEVEL_STYLES = dict(
            critical=dict(color="red", bold=True),
            error=dict(color="red"),
            warning=dict(color="yellow"),
            info=dict(color="green"),
            debug=dict(color="blue"),
        )
        coloredlogs.install(
            fmt=self.LOG_FORMAT,  # 使用你的自定义日志格式
            level=self.log_level,
            logger=logger,
            field_styles=dict(
                asctime=dict(color="yellow", bold=True),
                name=dict(color="blue", bold=True),
                filename=dict(color="cyan"),
                lineno=dict(color="green", bold=True),
                levelname=dict(color="magenta", bold=True),
            ),
            level_styles=dict(
                critical=dict(color="red", bold=True),
                error=dict(color="red"),
                warning=dict(color="yellow"),
                info=dict(color="green"),
                debug=dict(color="blue"),
            ),
        )
        logger.propagate = False
        return logger

    def set_names(self, cls_name: str | None = None, func_name: str | None = None):
        """
        Update the class and/or function context for the logger.

        Args:
            cls_name (str, optional): New class name.
            func_name (str, optional): New function name.
        """
        self.cls_name = cls_name
        self.func_name = func_name
        self.logger.name = self._get_logger_name()
        self.logger.propagate = False

    def set_level(self, log_level: str):
        """
        Update the logging level dynamically.

        Args:
            log_level (str): The new logging level.
        """
        self.log_level = self._get_log_level(log_level)
        self.logger.setLevel(self.log_level)
        coloredlogs.install(
            fmt=self.LOG_FORMAT,  # 使用你的自定义日志格式
            level=self.log_level,
            logger=self.logger,
            field_styles=dict(
                asctime=dict(color="yellow", bold=True),
                name=dict(color="blue", bold=True),
                filename=dict(color="cyan"),
                lineno=dict(color="green", bold=True),
                levelname=dict(color="magenta", bold=True),
            ),
            level_styles=dict(
                critical=dict(color="red", bold=True),
                error=dict(color="red"),
                warning=dict(color="yellow"),
                info=dict(color="green"),
                debug=dict(color="blue"),
            ),
            reconfigure=True,
        )
        self.logger.propagate = False

    def set_file(
        self,
        file: str,
        mode: str = "a",
        max_bytes: int = 10 * 1024 * 1024,
        backup_count: int = 5,
    ):
        """
        Configure file logging for the logger.

        Args:
            file (str): The name of the log file.
            mode (str, optional): File mode. Use 'a' for append or 'w' for overwrite. Defaults to "a".
            max_bytes (int, optional): Maximum size of a single log file before rotation. Defaults to 10MB.
            backup_count (int, optional): Number of backup files to keep. Defaults to 5.
        """
        if mode not in {"a", "w"}:
            raise BioatInvalidParameterError(
                f"Invalid file mode: {mode}. Choose 'a' (append) or 'w' (overwrite)."
            )
        # Ensure the directory exists
        directory = os.path.dirname(file)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
            self.logger.debug(f"Created directory for log file: {directory}")

        # Remove existing file handlers to prevent duplication
        self.logger.handlers = [
            h for h in self.logger.handlers if not isinstance(h, logging.FileHandler)
        ]

        # Use FileHandler for overwrite mode or RotatingFileHandler for append/rotation
        if mode == "w" or max_bytes <= 0:
            # Simple overwrite mode
            file_handler = logging.FileHandler(file, mode=mode)
        else:
            # Rotating mode
            file_handler = logging.handlers.RotatingFileHandler(
                file, mode=mode, maxBytes=max_bytes, backupCount=backup_count
            )

        # Set file log format
        file_handler.setFormatter(logging.Formatter(self.LOG_FORMAT))
        self.logger.addHandler(file_handler)
        self.logger.propagate = False

    def add_stream_handler(self):
        """
        Add a stream handler to the logger for writing logs to the console.
        """
        # 检查是否已有 StreamHandler
        if not any(isinstance(h, logging.StreamHandler) for h in self.logger.handlers):
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(logging.Formatter(self.LOG_FORMAT))
            self.logger.addHandler(stream_handler)
            self.logger.propagate = False


if __name__ == "__main__":
    # 初始化 LoggerManager
    lm = LoggerManager(log_level="debug")

    print("=== 基础日志测试 ===")
    lm.logger.debug("Debug message: Basic logging test.")
    lm.logger.info("Info message: Basic logging test.")
    lm.logger.warning("Warning message: Basic logging test.")
    lm.logger.error("Error message: Basic logging test.")
    lm.logger.critical("Critical message: Basic logging test.")

    print("\n=== 配置文件日志（追加模式） ===")
    lm.set_file("append_test.log", mode="a")
    lm.logger.info("Info message: Appending to 'append_test.log'.")
    lm.logger.warning("Warning message: Appending to 'append_test.log'.")
    lm.logger.error("Error message: Appending to 'append_test.log'.")

    print("\n=== 配置文件日志（覆盖模式） ===")
    lm.set_file("overwrite_test.log", mode="w")
    lm.logger.info("Info message: Overwriting 'overwrite_test.log'.")
    lm.logger.warning("Warning message: Overwriting 'overwrite_test.log'.")
    lm.logger.error("Error message: Overwriting 'overwrite_test.log'.")

    print("\n=== 配置文件日志（追加到不同文件） ===")
    # 追加另一个文件的同时不影响其他日志的写出
    lm.set_file("different_file.log", mode="a")
    lm.logger.info("Info message: Appending to 'different_file.log'.")
    lm.logger.warning("Warning message: Appending to 'different_file.log'.")
    lm.logger.error("Error message: Appending to 'different_file.log'.")

    print("\n=== 动态更新日志级别 ===")
    lm.set_level("WARNING")
    lm.logger.debug("Debug message: Should not appear (level=WARNING).")
    lm.logger.info("Info message: Should not appear (level=WARNING).")
    lm.logger.warning("Warning message: Should appear (level=WARNING).")
    lm.logger.error("Error message: Should appear (level=WARNING).")

    print("\n=== 动态更新类名和函数名 ===")
    lm.set_level("debug")
    lm.set_names(cls_name="TestClass", func_name="test_function")
    lm.logger.info("Info message: Logging with class and function names.")

    print("\n=== 多次调用文件日志写入 ===")
    lm.set_file("multi_file1.log", mode="w")
    lm.logger.info("First file: multi_file1.log.")
    lm.set_file("multi_file2.log", mode="w")
    lm.logger.info("Second file: multi_file2.log.")
    lm.set_file("multi_file3.log", mode="w")
    lm.logger.info("Third file: multi_file3.log.")

    print("\n=== 测试不合法的日志级别 ===")
    try:
        lm.set_level("fake_level")
    except BioatInvalidParameterError as err:
        print("Should raise FbtValueError")
        print(err)

    print("\n=== 测试重复调用 add_stream_handler ===")
    lm.add_stream_handler()  # 第一次调用
    lm.logger.info("First stream handler added.")
    lm.add_stream_handler()  # 第二次调用，不应重复添加
    lm.logger.info("No duplicate stream handler should be added.")

    print("\n=== 测试自动创建目录 ===")
    lm.set_file("logs/test_folder/append_test.log", mode="a")
    lm.logger.info(
        "This message should be written to 'append_test.log' inside 'logs/test_folder'."
    )
