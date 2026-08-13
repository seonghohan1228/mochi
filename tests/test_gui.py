import math

from mochi.gui import FULL_TURN_RAD, stop_angle_crossed


def test_stop_triggers_when_target_lies_inside_the_step() -> None:
    assert stop_angle_crossed(1.0, 1.2, 1.1)


def test_stop_ignores_target_beyond_the_step() -> None:
    assert not stop_angle_crossed(1.0, 1.2, 1.3)


def test_stop_ignores_target_just_passed() -> None:
    assert not stop_angle_crossed(1.0, 1.2, 0.99)


def test_stop_catches_target_across_the_wrap_at_full_turn() -> None:
    previous = FULL_TURN_RAD - 0.05
    advanced = previous + 0.1
    assert stop_angle_crossed(previous, advanced, 0.02)


def test_resuming_exactly_on_target_completes_a_full_revolution() -> None:
    target = 0.5 * math.pi
    assert not stop_angle_crossed(target, target + 0.1, target)
    assert stop_angle_crossed(target, target + FULL_TURN_RAD, target)


def test_zero_or_backward_step_never_stops() -> None:
    assert not stop_angle_crossed(1.0, 1.0, 1.0)
    assert not stop_angle_crossed(1.0, 0.9, 0.95)


def test_step_longer_than_one_revolution_always_stops() -> None:
    assert stop_angle_crossed(0.0, 1.5 * FULL_TURN_RAD, 4.0)
