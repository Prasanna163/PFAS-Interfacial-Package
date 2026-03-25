import os
import sys


def _clear_screen():
    os.system("cls" if os.name == "nt" else "clear")


def _read_key_windows():
    import msvcrt

    key = msvcrt.getwch()
    if key in ("\x00", "\xe0"):
        key = msvcrt.getwch()
        mapping = {
            "H": "up",
            "P": "down",
            "K": "left",
            "M": "right",
        }
        return mapping.get(key, "other")
    if key == "\r":
        return "enter"
    if key == "\x03":
        raise KeyboardInterrupt
    return key


def _read_key_unix():
    import termios
    import tty

    stream = sys.stdin
    fileno = stream.fileno()
    original = termios.tcgetattr(fileno)
    try:
        tty.setraw(fileno)
        first = stream.read(1)
        if first == "\x03":
            raise KeyboardInterrupt
        if first in ("\r", "\n"):
            return "enter"
        if first == "\x1b":
            second = stream.read(1)
            third = stream.read(1)
            if second == "[":
                return {
                    "A": "up",
                    "B": "down",
                    "C": "right",
                    "D": "left",
                }.get(third, "other")
            return "other"
        return first
    finally:
        termios.tcsetattr(fileno, termios.TCSADRAIN, original)


def read_key():
    return _read_key_windows() if os.name == "nt" else _read_key_unix()


def choose_option(title, options, footer=None):
    if not options:
        raise ValueError("choose_option requires at least one option.")

    selected = 0
    while True:
        _clear_screen()
        print(title)
        print("-" * len(title))
        for index, option in enumerate(options):
            marker = ">" if index == selected else " "
            print(f"{marker} {option}")
        print("")
        if footer:
            print(footer)
        print("Use Up/Down arrows and Enter to select.")

        key = read_key()
        if key == "up":
            selected = (selected - 1) % len(options)
        elif key == "down":
            selected = (selected + 1) % len(options)
        elif key == "enter":
            _clear_screen()
            return selected


def prompt_text(label, default=None, allow_empty=False):
    while True:
        suffix = f" [{default}]" if default is not None else ""
        value = input(f"{label}{suffix}: ").strip()
        if not value and default is not None:
            value = str(default)
        if value or allow_empty:
            return value
        print("A value is required.")


def prompt_float(label, default=None, min_value=None):
    while True:
        raw = prompt_text(label, default=default)
        try:
            value = float(raw)
        except ValueError:
            print("Please enter a valid number.")
            continue
        if min_value is not None and value < min_value:
            print(f"Value must be >= {min_value}.")
            continue
        return value


def prompt_int(label, default=None):
    while True:
        raw = prompt_text(label, default=default)
        try:
            return int(raw)
        except ValueError:
            print("Please enter a valid integer.")


def ensure_interactive_terminal():
    if not sys.stdin.isatty() or not sys.stdout.isatty():
        raise RuntimeError("Interactive mode requires a terminal (TTY).")


__all__ = [
    "choose_option",
    "ensure_interactive_terminal",
    "prompt_float",
    "prompt_int",
    "prompt_text",
]
