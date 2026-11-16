import sys
if len(sys.argv)>1:
    from .cli import main
    main(sys.argv[1:])
else:
    try:
        from .gui import launch
    except ImportError as exc:
        raise SystemExit(f"GUI dependency is missing: {exc}\nInstall with: pip install -r requirements.txt")
    raise SystemExit(launch())
