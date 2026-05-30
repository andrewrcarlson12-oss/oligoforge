"""OligoForge desktop launcher.

Frozen into a double-click executable by PyInstaller (see oligoforge.spec):
it boots the local server, waits until it is listening, then opens the cockpit
in a native window if pywebview is installed, otherwise in the default browser.

Run from source to test the same flow:  python launcher.py
"""
import contextlib
import socket
import threading
import time
import webbrowser

HOST = "127.0.0.1"
APP_VERSION = "v1.1.1"                         # bumped on each change; matches the release tag
GITHUB_REPO = "andrewrcarlson12-oss/oligoforge"  # <- set to your repo for update checks


def _ver(v):
    return tuple(int(x) for x in v.lstrip("vV").split(".") if x.isdigit())


def check_for_update():
    """Quietly tell the user if a newer GitHub release exists. Never blocks or fails loudly."""
    try:
        import json
        import urllib.request
        url = "https://api.github.com/repos/%s/releases/latest" % GITHUB_REPO
        req = urllib.request.Request(url, headers={"Accept": "application/vnd.github+json",
                                                   "User-Agent": "OligoForge"})
        with urllib.request.urlopen(req, timeout=4) as r:
            data = json.load(r)
        latest = data.get("tag_name") or ""
        if latest and _ver(latest) > _ver(APP_VERSION):
            link = data.get("html_url") or ("https://github.com/%s/releases/latest" % GITHUB_REPO)
            print("\n  *** Update available: %s (you have %s) ***" % (latest, APP_VERSION))
            print("  Download the new one-click app: %s\n" % link)
    except Exception:
        pass


def _free_port(preferred=8111):
    """Return the preferred port if nothing is listening on it, else the next free one."""
    for p in [preferred] + list(range(8112, 8140)):
        with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
            if s.connect_ex((HOST, p)) != 0:   # connect failed -> nobody is listening -> free
                return p
    return preferred


def _wait_until_up(port, timeout=30.0):
    end = time.time() + timeout
    while time.time() < end:
        with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
            if s.connect_ex((HOST, port)) == 0:
                return True
        time.sleep(0.15)
    return False


def _serve(port):
    import uvicorn
    from app import app
    uvicorn.run(app, host=HOST, port=port, log_level="warning")


def main():
    threading.Thread(target=check_for_update, daemon=True).start()
    port = _free_port()
    url = "http://localhost:%d" % port  # localhost (Firefox dislikes the bare IP); server still binds 127.0.0.1

    threading.Thread(target=_serve, args=(port,), daemon=True).start()
    if not _wait_until_up(port):
        print("OligoForge could not start its local server. "
              "Make sure nothing else is using the port, then try again.")
        try:
            input("Press Enter to close this window... ")
        except EOFError:
            pass
        return

    # Best experience: a native window. Falls back to the browser if pywebview
    # (or its system web runtime) isn't available.
    try:
        import webview
        webview.create_window("OligoForge", url, width=1180, height=900)
        webview.start()
        return  # window closed -> process exits -> daemon server stops
    except Exception:
        pass

    print("\n  OligoForge is running.")
    print("  Opening %s in your browser." % url)
    print("  Keep this window open while you work. Close it to stop OligoForge.\n")
    with contextlib.suppress(Exception):
        webbrowser.open(url)
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
