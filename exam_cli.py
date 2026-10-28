"""Chạy một bài JSON của PPS trong terminal, không cần khởi tạo cửa sổ Tk."""
import json
import sys
from pathlib import Path

from pps.gui import dispatch


def main(argv=None):
    argv=sys.argv[1:] if argv is None else argv
    if len(argv)!=1:
        print('Cách dùng: python exam_cli.py duong_dan_bai.json',file=sys.stderr)
        return 2
    path=Path(argv[0])
    try:
        task=json.loads(path.read_text(encoding='utf-8-sig'))
        report=dispatch(task['tool'],task['method'],task.get('values',{}),
                        task.get('x',''),task.get('y',''))
        print(report.text())
        return 0
    except (OSError,UnicodeError,ValueError,KeyError,TypeError) as exc:
        print(f'Không tính được: {exc}',file=sys.stderr)
        return 1


if __name__=='__main__':raise SystemExit(main())
