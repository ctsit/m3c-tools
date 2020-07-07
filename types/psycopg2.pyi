from typing import List, Optional, Tuple, TextIO


class ConnectionFactory:
    ...


class CursorFactory:
    ...


class Connection:
    def __enter__(self) -> None:
        ...

    def __exit__(self, type, value, traceback) -> None:
        ...

    def close(self) -> None:
        ...

    def cursor(self) -> "Cursor":
        ...


class Cursor:
    rowcount: int

    def __enter__(self) -> "Cursor":
        ...

    def __exit__(self, type, value, traceback) -> None:
        ...

    def __iter__(self):
        ...

    def copy_from(self, src: TextIO, table: str,
                  columns: Optional[Tuple] = None) -> None:
        ...

    def execute(self, query: str, params: Optional[Tuple] = None) -> None:
        ...

    def fetchone(self) -> List:
        ...


def connect(
    dsn: Optional[str] = None,
    connection_factory: Optional[ConnectionFactory] = None,
    cursor_factory: Optional[CursorFactory] = None,
    **kwargs
) -> Connection:
    ...
