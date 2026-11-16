from __future__ import annotations

import math
from typing import Iterable

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPen, QBrush, QFont, QPolygonF
from PySide6.QtWidgets import (
    QGraphicsEllipseItem,
    QGraphicsItem,
    QGraphicsLineItem,
    QGraphicsPolygonItem,
    QGraphicsScene,
    QGraphicsSimpleTextItem,
    QGraphicsView,
)


class _NodeItem(QGraphicsEllipseItem):
    def __init__(self, view: "RVESketchView", node_id: int, x: float, y: float):
        r = 2.4
        super().__init__(-r, -r, 2 * r, 2 * r)
        self.view = view
        self.node_id = node_id
        self.setPos(x, -y)
        self.setBrush(QBrush(QColor("#FFFFFF")))
        self.setPen(QPen(QColor("#2F80ED"), 1.15))
        self.setFlags(
            QGraphicsItem.GraphicsItemFlag.ItemIsMovable
            | QGraphicsItem.GraphicsItemFlag.ItemIsSelectable
            | QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges
        )
        self.setZValue(20)
        self.label = QGraphicsSimpleTextItem(f"N{node_id}")
        self.label.setFont(QFont("Segoe UI", 8))
        self.label.setBrush(QBrush(QColor("#172033")))
        self.label.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
        self.label.setZValue(25)
        self.view._scene.addItem(self.label)
        self._place_label()

    def _place_label(self):
        p = self.pos()
        self.label.setPos(p.x() + 4.0, p.y() - 14.0)

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionChange and self.view.snap_enabled and not self.view._loading:
            p = value
            g = self.view.grid_spacing
            if g > 0:
                value = QPointF(round(p.x() / g) * g, round(p.y() / g) * g)
        out = super().itemChange(change, value)
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            self._place_label()
            self.view._node_moved(self)
        if change == QGraphicsItem.GraphicsItemChange.ItemSelectedHasChanged:
            selected = bool(value)
            self.setBrush(QBrush(QColor("#DCEBFF" if selected else "#FFFFFF")))
        return out

    def remove_label(self):
        if self.label.scene() is not None:
            self.view._scene.removeItem(self.label)


class _WallItem(QGraphicsLineItem):
    def __init__(self, view: "RVESketchView", wall_id: int, node_i: int, node_j: int, weight: float = 1.0, stack_id: int = 1):
        super().__init__()
        self.view = view
        self.wall_id = wall_id
        self.node_i = node_i
        self.node_j = node_j
        self.weight = float(weight)
        self.stack_id = int(stack_id)
        self.setPen(QPen(QColor("#172033"), 1.0))
        self.setFlags(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setZValue(10)
        self.label = QGraphicsSimpleTextItem(f"W{wall_id}")
        self.label.setFont(QFont("Segoe UI", 8))
        self.label.setBrush(QBrush(QColor("#6B7280")))
        self.label.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
        self.label.setZValue(15)
        self.view._scene.addItem(self.label)
        # Arrowhead at mid-wall pointing from Node I to Node J: the positive local 1̄ axis.
        # A child item (removed with the wall), constant size on screen, never takes the mouse.
        self.arrow = QGraphicsPolygonItem(QPolygonF([QPointF(6.0, 0.0), QPointF(-4.5, -4.2), QPointF(-4.5, 4.2)]), self)
        self.arrow.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations, True)
        self.arrow.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
        self.arrow.setPen(QPen(Qt.PenStyle.NoPen))
        self.arrow.setBrush(QBrush(QColor("#172033")))
        self.update_geometry()

    def update_geometry(self):
        a = self.view.nodes[self.node_i - 1].pos()
        b = self.view.nodes[self.node_j - 1].pos()
        self.setLine(a.x(), a.y(), b.x(), b.y())
        self.arrow.setPos((a.x() + b.x()) / 2, (a.y() + b.y()) / 2)
        self.arrow.setRotation(math.degrees(math.atan2(b.y() - a.y(), b.x() - a.x())))
        self.arrow.setVisible(self.view.show_wall_direction and (a - b).manhattanLength() > 1e-12)
        self.setToolTip(f"W{self.wall_id}: N{self.node_i} → N{self.node_j}\nThe local 1̄ axis points from "
                        f"Node I (N{self.node_i}) to Node J (N{self.node_j}).")
        self.label.setPos((a.x() + b.x()) / 2 + 4, (a.y() + b.y()) / 2 + 4)
        if self.view.show_dimensions:
            dx=b.x()-a.x(); dy_model=-(b.y()-a.y())
            length=math.hypot(dx,dy_model); angle=math.degrees(math.atan2(dy_model,dx))
            self.label.setText(f"W{self.wall_id}  {length:.4g} mm  {angle:.3g}°")
        else:
            self.label.setText(f"W{self.wall_id}")

    def itemChange(self, change, value):
        out=super().itemChange(change,value)
        if change == QGraphicsItem.GraphicsItemChange.ItemSelectedHasChanged:
            selected=bool(value)
            self.setPen(QPen(QColor("#2F80ED" if selected else "#172033"), 1.6 if selected else 1.0))
            self.arrow.setBrush(QBrush(QColor("#2F80ED" if selected else "#172033")))
        return out

    def remove_label(self):
        if self.label.scene() is not None:
            self.view._scene.removeItem(self.label)


class RVESketchView(QGraphicsView):
    """Small CAD-like view for defining the 2-D RVE graph.

    Scene Y is inverted relative to model Y, so the user sees a conventional Cartesian
    x1-x2 plane while all exported coordinates remain in the solver convention.
    """

    geometryEdited = Signal()
    latticeVectorDefined = Signal(str, float, float)
    statusMessage = Signal(str)

    MODES = {"select", "node", "wall", "delete", "a1", "a2"}

    def __init__(self, parent=None):
        # QGraphicsView does not take ownership of a scene passed to its constructor.
        # Keep the scene as a persistent, view-owned QObject so Python cannot garbage-
        # collect it after __init__ returns (which would make self.scene() return None).
        super().__init__(parent)
        self._scene = QGraphicsScene(self)
        self.setScene(self._scene)
        self.nodes: list[_NodeItem] = []
        self.walls: list[_WallItem] = []
        self.mode = "select"
        self.grid_spacing = 5.0
        self.snap_enabled = True
        self._pending_node: int | None = None
        self._loading = False
        self.show_dimensions = False
        self.show_wall_direction = True
        self._vector_items = []
        self.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        self.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self.setBackgroundBrush(QBrush(QColor("#FCFDFE")))
        self.setSceneRect(QRectF(-250, -250, 500, 500))

    def set_mode(self, mode: str):
        if mode not in self.MODES:
            raise ValueError(f"Unknown sketch mode {mode}")
        self.mode = mode
        self._pending_node = None
        self._clear_pending_highlight()
        self.setDragMode(QGraphicsView.DragMode.RubberBandDrag if mode == "select" else QGraphicsView.DragMode.NoDrag)
        messages = {
            "select": "Select/drag nodes. Geometry values update automatically.",
            "node": "Click the grid to place a node.",
            "wall": "Click two existing nodes to draw a wall.",
            "delete": "Click a node or wall to delete it.",
            "a1": "Click two nodes: the vector from the first to the second becomes a1.",
            "a2": "Click two nodes: the vector from the first to the second becomes a2.",
        }
        self.statusMessage.emit(messages[mode])

    def set_grid_spacing(self, spacing: float):
        self.grid_spacing = max(float(spacing), 1e-9)
        self.viewport().update()

    def set_snap_enabled(self, enabled: bool):
        self.snap_enabled = bool(enabled)

    def drawBackground(self, painter: QPainter, rect: QRectF):
        super().drawBackground(painter, rect)
        g = self.grid_spacing
        if g <= 0:
            return
        left = math.floor(rect.left() / g) * g
        top = math.floor(rect.top() / g) * g
        lines = []
        x = left
        while x <= rect.right():
            lines.append((QPointF(x, rect.top()), QPointF(x, rect.bottom())))
            x += g
        y = top
        while y <= rect.bottom():
            lines.append((QPointF(rect.left(), y), QPointF(rect.right(), y)))
            y += g
        painter.setPen(QPen(QColor("#E9EEF5"), 0.0))
        for a, b in lines:
            painter.drawLine(a, b)
        painter.setPen(QPen(QColor("#B8C5D6"), 0.0))
        painter.drawLine(QPointF(rect.left(), 0), QPointF(rect.right(), 0))
        painter.drawLine(QPointF(0, rect.top()), QPointF(0, rect.bottom()))

    def wheelEvent(self, event):
        factor = 1.18 if event.angleDelta().y() > 0 else 1 / 1.18
        self.scale(factor, factor)

    def mousePressEvent(self, event):
        if event.button() != Qt.MouseButton.LeftButton or self.mode == "select":
            super().mousePressEvent(event)
            return
        p = self.mapToScene(event.position().toPoint())
        if self.mode == "node":
            x, y = self._snapped_model_point(p)
            self.add_node(x, y)
            self.geometryEdited.emit()
            return
        if self.mode in ("wall", "a1", "a2"):
            node_id = self._nearest_node_id(p)
            if node_id is None:
                self.statusMessage.emit("Click close to an existing node.")
                return
            self._two_node_action(node_id)
            return
        if self.mode == "delete":
            item = self.itemAt(event.position().toPoint())
            if item is not None and isinstance(item.parentItem(), _WallItem):
                item = item.parentItem()  # the wall's direction arrowhead
            if isinstance(item, _NodeItem):
                self.delete_node(item.node_id)
                self.geometryEdited.emit()
            elif isinstance(item, _WallItem):
                self.delete_wall(item.wall_id)
                self.geometryEdited.emit()
            return

    def _two_node_action(self, node_id: int):
        if self._pending_node is None:
            self._pending_node = node_id
            self.nodes[node_id - 1].setBrush(QBrush(QColor("#FDE6B7")))
            self.statusMessage.emit(f"First node N{node_id} selected; choose the second node.")
            return
        first = self._pending_node
        self._pending_node = None
        self._clear_pending_highlight()
        if first == node_id:
            self.statusMessage.emit("Choose two different nodes.")
            return
        if self.mode == "wall":
            if self._wall_exists(first, node_id):
                self.statusMessage.emit("That wall already exists.")
                return
            self.add_wall(first, node_id)
            self.geometryEdited.emit()
        else:
            a = self.model_position(first)
            b = self.model_position(node_id)
            dx, dy = b[0] - a[0], b[1] - a[1]
            self.latticeVectorDefined.emit(self.mode, dx, dy)
            self.statusMessage.emit(f"{self.mode} = [{dx:.6g}, {dy:.6g}]")

    def _clear_pending_highlight(self):
        for n in self.nodes:
            if not n.isSelected():
                n.setBrush(QBrush(QColor("#FFFFFF")))

    def _snapped_model_point(self, p: QPointF) -> tuple[float, float]:
        x, y = p.x(), -p.y()
        if self.snap_enabled:
            g = self.grid_spacing
            x, y = round(x / g) * g, round(y / g) * g
        return float(x), float(y)

    def _nearest_node_id(self, p: QPointF) -> int | None:
        if not self.nodes:
            return None
        # Map a ~12 pixel pick radius to scene units so picking remains comfortable after zoom.
        p2 = self.mapToScene(self.mapFromScene(p).x() + 12, self.mapFromScene(p).y())
        radius = max(abs(p2.x() - p.x()), 3.0)
        best = None
        best_d = float("inf")
        for n in self.nodes:
            d = math.hypot(n.pos().x() - p.x(), n.pos().y() - p.y())
            if d < best_d and d <= radius:
                best, best_d = n.node_id, d
        return best

    def _wall_exists(self, ni: int, nj: int) -> bool:
        key = {ni, nj}
        return any({w.node_i, w.node_j} == key for w in self.walls)

    def _node_moved(self, node: _NodeItem):
        if self._loading:
            return
        for w in self.walls:
            if w.node_i == node.node_id or w.node_j == node.node_id:
                w.update_geometry()
        self.geometryEdited.emit()

    def model_position(self, node_id: int) -> tuple[float, float]:
        p = self.nodes[node_id - 1].pos()
        return float(p.x()), float(-p.y())

    def add_node(self, x: float, y: float):
        item = _NodeItem(self, len(self.nodes) + 1, float(x), float(y))
        self._scene.addItem(item)
        self.nodes.append(item)
        self._expand_scene_rect()

    def add_wall(self, node_i: int, node_j: int, weight: float = 1.0, stack_id: int = 1):
        if node_i < 1 or node_j < 1 or node_i > len(self.nodes) or node_j > len(self.nodes):
            return
        item = _WallItem(self, len(self.walls) + 1, int(node_i), int(node_j), weight, stack_id)
        self._scene.addItem(item)
        self.walls.append(item)

    def delete_wall(self, wall_id: int):
        if wall_id < 1 or wall_id > len(self.walls):
            return
        w = self.walls.pop(wall_id - 1)
        w.remove_label()
        self._scene.removeItem(w)
        for i, wall in enumerate(self.walls, 1):
            wall.wall_id = i
            wall.label.setText(f"W{i}")

    def delete_node(self, node_id: int):
        if node_id < 1 or node_id > len(self.nodes):
            return
        # Capture data, remove target, remap graph, and rebuild to guarantee contiguous 1-based IDs.
        old_nodes = [self.model_position(i + 1) for i in range(len(self.nodes))]
        old_walls = [(w.node_i, w.node_j, w.weight, w.stack_id) for w in self.walls]
        del old_nodes[node_id - 1]
        new_walls = []
        for ni, nj, wt, sid in old_walls:
            if node_id in (ni, nj):
                continue
            ni2 = ni - 1 if ni > node_id else ni
            nj2 = nj - 1 if nj > node_id else nj
            new_walls.append((ni2, nj2, wt, sid))
        self.set_geometry(old_nodes, new_walls, fit=False)

    def clear_geometry(self):
        self._loading = True
        try:
            for w in self.walls:
                w.remove_label()
            for n in self.nodes:
                n.remove_label()
            self._scene.clear()
            self.nodes = []
            self.walls = []
            self._vector_items = []
            self._pending_node = None
        finally:
            self._loading = False

    def set_geometry(self, nodes: Iterable[tuple[float, float]], walls: Iterable[tuple], fit: bool = True):
        self.clear_geometry()
        self._loading = True
        try:
            for x, y in nodes:
                self.add_node(float(x), float(y))
            for row in walls:
                ni, nj = int(row[0]), int(row[1])
                wt = float(row[2]) if len(row) >= 3 else 1.0
                sid = int(row[3]) if len(row) >= 4 else 1
                self.add_wall(ni, nj, wt, sid)
        finally:
            self._loading = False
        self._expand_scene_rect()
        if fit:
            self.fit_geometry()

    def export_geometry(self):
        nodes = [self.model_position(i + 1) for i in range(len(self.nodes))]
        walls = [(w.node_i, w.node_j, w.weight, w.stack_id) for w in self.walls]
        return nodes, walls

    def set_lattice_vectors(self, a1, a2):
        """Draw a1/a2 and the translated parallel edges of the periodic cell.

        The extra two lines make the periodic parallelogram immediately visible:
        a line parallel to a1 is drawn from a2 to a1+a2, and a line parallel to
        a2 is drawn from a1 to a1+a2.
        """
        for item in self._vector_items:
            try:
                if item.scene() is not None:
                    self._scene.removeItem(item)
            except RuntimeError:
                pass
        self._vector_items = []
        try:
            u = (float(a1[0]), float(a1[1]))
            v = (float(a2[0]), float(a2[1]))
        except Exception:
            self._expand_scene_rect()
            return
        valid_u = all(math.isfinite(x) for x in u) and math.hypot(*u) > 1e-12
        valid_v = all(math.isfinite(x) for x in v) and math.hypot(*v) > 1e-12

        def add_line(x1, y1, x2, y2, color, style=Qt.PenStyle.DashLine, width=1.25):
            line = QGraphicsLineItem(x1, -y1, x2, -y2)
            line.setPen(QPen(QColor(color), width, style))
            line.setZValue(4)
            self._scene.addItem(line)
            self._vector_items.append(line)
            return line

        def add_label(text, x, y, color):
            label = QGraphicsSimpleTextItem(text)
            label.setFont(QFont("Segoe UI", 8))
            label.setBrush(QBrush(QColor(color)))
            label.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
            label.setPos(x + 4, -y + 4)
            label.setZValue(5)
            self._scene.addItem(label)
            self._vector_items.append(label)

        if valid_u:
            add_line(0.0, 0.0, u[0], u[1], "#2F80ED")
            add_label("a1", u[0], u[1], "#2F80ED")
        if valid_v:
            add_line(0.0, 0.0, v[0], v[1], "#7A5AF8")
            add_label("a2", v[0], v[1], "#7A5AF8")
        if valid_u and valid_v:
            # Opposite periodic-cell edges, explicitly parallel to a1 and a2.
            p3 = (u[0] + v[0], u[1] + v[1])
            add_line(v[0], v[1], p3[0], p3[1], "#2F80ED", Qt.PenStyle.DotLine, 1.0)
            add_line(u[0], u[1], p3[0], p3[1], "#7A5AF8", Qt.PenStyle.DotLine, 1.0)
        self._expand_scene_rect()


    def set_show_dimensions(self, enabled: bool):
        self.show_dimensions=bool(enabled)
        for wall in self.walls:
            wall.update_geometry()

    def set_show_wall_direction(self, enabled: bool):
        """Arrowheads Node I -> Node J (positive local 1̄ axis) on every wall."""
        self.show_wall_direction=bool(enabled)
        for wall in self.walls:
            wall.update_geometry()

    def selected_wall_ids(self) -> list[int]:
        return [w.wall_id for w in self.walls if w.isSelected()]

    def highlight_walls(self, wall_ids: Iterable[int]):
        wanted={int(i) for i in wall_ids}
        for w in self.walls:
            w.setSelected(w.wall_id in wanted)
        if wanted:
            self.statusMessage.emit("Highlighted " + ", ".join(f"W{i}" for i in sorted(wanted)))

    def constrain_wall(self, wall_id: int, kind: str, value: float | None = None) -> bool:
        """Apply a simple sketch constraint to a wall by moving its second node.

        These are editing aids, not solver constraints: the resulting coordinates are
        written back to the normal node/wall representation used by the solver.
        """
        if wall_id < 1 or wall_id > len(self.walls):
            return False
        w=self.walls[wall_id-1]
        a=self.nodes[w.node_i-1].pos(); b=self.nodes[w.node_j-1].pos()
        dx=b.x()-a.x(); dy=b.y()-a.y()
        if kind == "horizontal":
            b2=QPointF(b.x(),a.y())
        elif kind == "vertical":
            b2=QPointF(a.x(),b.y())
        elif kind == "length":
            target=float(value or 0.0)
            L=math.hypot(dx,dy)
            if target <= 0 or L <= 1e-12:
                return False
            b2=QPointF(a.x()+dx*target/L,a.y()+dy*target/L)
        else:
            return False
        old=self._loading; self._loading=True
        try:
            self.nodes[w.node_j-1].setPos(b2)
            self.nodes[w.node_j-1]._place_label()
            for wall in self.walls:
                if wall.node_i == w.node_j or wall.node_j == w.node_j:
                    wall.update_geometry()
        finally:
            self._loading=old
        self.geometryEdited.emit()
        return True

    def fit_geometry(self):
        if not self.nodes:
            self.resetTransform()
            self.setSceneRect(QRectF(-100, -100, 200, 200))
            return
        rect = self._scene.itemsBoundingRect().adjusted(-25, -25, 25, 25)
        if rect.width() < 1:
            rect.adjust(-50, 0, 50, 0)
        if rect.height() < 1:
            rect.adjust(0, -50, 0, 50)
        self.setSceneRect(rect)
        self.fitInView(rect, Qt.AspectRatioMode.KeepAspectRatio)

    def _expand_scene_rect(self):
        if not self.nodes:
            return
        rect = self._scene.itemsBoundingRect().adjusted(-100, -100, 100, 100)
        self.setSceneRect(self.sceneRect().united(rect))
