import logging
import logging.handlers
import os
import threading

import coloredlogs
import json_log_formatter

from bioat.exceptions import BioatInvalidParameterError

__all__ = ["LoggerManager"]

_logger_cache = {}
_logger_lock = threading.Lock()

class LoggerManager:
    LOG_FORMAT = "%(asctime)s.%(msecs)03d - [%(name)s] - %(filename)s[line:%(lineno)4d] - %(levelname)+8s: %(message)s"
    DEFAULT_LEVEL = os.getenv("BIOAT_LOG_LEVEL", "INFO").upper()

    def __init__(
        self,
        log_level: str = DEFAULT_LEVEL,
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
        self.logger = self._get_or_create_logger()
        # self.logger.propagate = False # TODO

    @staticmethod
    def get_logger(name: str, level: str = DEFAULT_LEVEL) -> logging.Logger:
        """直接返回 logger 实例，用于标准 logging 接口"""
        return LoggerManager(mod_name=name, log_level=level).logger

    def _get_log_level(self, log_level: str) -> int:
        if log_level.upper() not in logging._nameToLevel:
            raise BioatInvalidParameterError(
                f"Invalid log level: {log_level}. "
                "Choose from CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET."
            )
        return logging._nameToLevel[log_level.upper()]

    def _get_logger_name(self) -> str:
        parts = [self.mod_name]
        if self.cls_name:
            parts.append(self.cls_name)
        if self.func_name:
            parts.append(self.func_name)
        return ".".join(parts)

    def _get_or_create_logger(self) -> logging.Logger:
        name = self._get_logger_name()
        with _logger_lock:
            if name in _logger_cache:
                return _logger_cache[name]

            logger = logging.getLogger(name)
            logger.setLevel(self.log_level)
            logger.propagate = False

            if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
                coloredlogs.install(
                    fmt=self.LOG_FORMAT,
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
            _logger_cache[name] = logger
            return logger

    def set_names(self, cls_name: str | None = None, func_name: str | None = None):
        """重新设置类名和函数名，会重新绑定 logger"""
        self.cls_name = cls_name
        self.func_name = func_name
        # 会根据新的 name 自动获取缓存或新建 logger（线程安全）
        # 可能只是应用了同名 logger 的配置,并不会实例化新的 logger而是指向同名旧 logger
        # 也可能 logger 会被新的 logger 实例替换
        self.logger = self._get_or_create_logger()

    def set_level(self, log_level: str):
        """更新日志等级"""
        self.log_level = self._get_log_level(log_level)
        self.logger.setLevel(self.log_level)

    def mute(self):
        """将日志等级设置为 NOTSET"""
        self.set_level(log_level="NOTSET")  # NOTSET is higher than CRITICAL

    def set_file(
        self,
        file: str,
        mode: str = "a",
        max_bytes: int = 10 * 1024 * 1024,
        backup_count: int = 5,
    ):
        """设置普通文件日志"""
        if mode not in {"a", "w"}:
            raise BioatInvalidParameterError("Mode must be 'a' or 'w'")
        os.makedirs(os.path.dirname(file), exist_ok=True)

        self._remove_handler_type(logging.FileHandler)

        if mode == "w" or max_bytes <= 0:
            handler = logging.FileHandler(file, mode=mode)
        else:
            handler = logging.handlers.RotatingFileHandler(
                file, mode=mode, maxBytes=max_bytes, backupCount=backup_count
            )

        handler.setFormatter(logging.Formatter(self.LOG_FORMAT))
        self.logger.addHandler(handler)

    def set_json_file(self, file: str, mode="a"):
        """设置 JSON 格式的文件日志输出"""
        os.makedirs(os.path.dirname(file), exist_ok=True)
        self._remove_handler_type(logging.FileHandler)

        handler = logging.FileHandler(file, mode=mode)
        formatter = json_log_formatter.JSONFormatter()
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

    def _remove_handler_type(self, handler_type):
        """内部方法：移除已有的特定 handler 类型"""
        self.logger.handlers = [
            h for h in self.logger.handlers if not isinstance(h, handler_type)
        ]

    def add_stream_handler(self):
        """添加 console handler（不会重复添加）"""
        if not any(isinstance(h, logging.StreamHandler) for h in self.logger.handlers):
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(logging.Formatter(self.LOG_FORMAT))
            self.logger.addHandler(stream_handler)


if __name__ == "__main__":
    import os
    import tempfile

    print("=== LoggerManager 内部测试开始 ===")

    # 测试 logger 创建与等级
    lm = LoggerManager(mod_name="test.module", log_level="DEBUG")
    assert isinstance(lm.logger, logging.Logger)
    assert lm.logger.level == logging.DEBUG
    print("✅ logger 创建与等级 OK")

    # 测试 set_names
    lm.set_names(cls_name="TestClass", func_name="test_func")
    assert "TestClass" in lm.logger.name
    assert "test_func" in lm.logger.name
    print("✅ set_names OK:", lm.logger.name)

    # 测试 set_level
    lm.set_level("INFO")
    assert lm.logger.level == logging.INFO
    lm.logger.info("✅This is an info message.")
    lm.logger.debug("❌This debug message should not appear.")
    lm.logger.warning("✅This is a warning message.")
    lm.logger.error("✅This is an error message.")
    lm.logger.critical("✅This is a critical message.")
    lm.set_level("DEBUG")
    assert lm.logger.level == logging.DEBUG
    lm.logger.debug("✅This debug message should now appear.")
    lm.logger.info("✅This is another info message after level change.")
    print("✅ set_level OK")

    # 测试非法等级
    try:
        lm.set_level("FAKE_LEVEL")
    except BioatInvalidParameterError as e:
        print("✅ 错误等级捕获 OK:", type(e), e)
    else:
        raise AssertionError("❌ 错误等级未正确捕获")

    # 测试文件日志
    with tempfile.TemporaryDirectory() as tmpdir:
        file_path = os.path.join(tmpdir, "test.log")
        lm.set_file(file_path, mode="w")
        lm.logger.warning("test file log entry")
        with open(file_path) as f:
            content = f.read()
            assert "test file log entry" in content
    print("✅ 文件日志写入 OK:", file_path)

    # 测试 JSON 日志
    with tempfile.TemporaryDirectory() as tmpdir:
        json_path = os.path.join(tmpdir, "log.json")
        lm.set_json_file(json_path)
        lm.logger.error("json log message")
        with open(json_path) as f:
            line = f.readline()
            assert '"message": "json log message"' in line
    print("✅ JSON 文件日志 OK:", json_path)

    # 测试 add_stream_handler 不重复
    # 仅统计 StreamHandler 类型的 handler 数量
    def count_stream_handlers(logger):
        return sum(1 for h in logger.handlers if isinstance(h, logging.StreamHandler))

    count_before = count_stream_handlers(lm.logger)
    lm.add_stream_handler()
    lm.add_stream_handler()
    count_after = count_stream_handlers(lm.logger)

    # 至多只能添加一个额外的 StreamHandler（如果 coloredlogs 没装）
    assert count_after == count_before, "Unexpected number of StreamHandlers"
    print(
        "✅ StreamHandler 不重复添加 OK:",
        f"StreamHandlers: before={count_before}, after={count_after}",
    )

    # 静态方式测试
    static_logger = LoggerManager.get_logger("static.test", level="WARNING")
    assert isinstance(static_logger, logging.Logger)
    assert "static.test" in static_logger.name
    print("✅ get_logger 静态方法 OK:", static_logger.name)

    lm = LoggerManager(mod_name="bioat.lib.align", log_level="DEBUG")
    lm.logger.info("Hello Logger!")

    # 修改类名函数名
    lm.set_names(cls_name="MyClass", func_name="run")
    lm.logger.debug("With class+func name.")

    # 使用 JSON 文件日志
    lm.set_json_file("logs/bioat_log.json")

    # 静态方式快速获取 logger
    logger = LoggerManager.get_logger("bioat.io", "INFO")
    logger.info("Hello from static logger")
    # 初始化 LoggerManager
    lm = LoggerManager(log_level="debug")
    print("🎉 所有 LoggerManager 内部测试通过！")
    print("_logger_cache:", _logger_cache)  # 打印缓存的 logger 名称
