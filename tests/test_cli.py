from click.testing import CliRunner

from oelsmann24_verticallandmotion.cli import main


def test_cli_help() -> None:
    """The installed Click entry point exposes its documented options."""
    result = CliRunner().invoke(main, ["--help"])

    assert result.exit_code == 0
    assert "Usage: main [OPTIONS]" in result.output
    assert "--input-vlm-reconstruction-file" in result.output
    assert "--output-lslr-file" in result.output
    assert "--vlm-rate-mode" in result.output
    assert "converge-to-gia" in result.output
    assert "constant" in result.output
    assert "--gia-only" in result.output
