#![feature(cmp_minmax)]

use eframe::egui::{self, Color32, FontId};
use fmtastic::*;
use lute::atom_info::ATOM_DATA;
use lute::disp::svg::svg_atom_color;
use lute::prelude::*;
use petgraph::visit::*;
use std::fmt::Write;
use tracing::level_filters::LevelFilter;
use tracing_subscriber::{filter, prelude::*, Registry};

#[cfg(not(target_family = "wasm"))]
fn main() -> eframe::Result {
    let mut arena = Arena::<u32>::new();
    let mut current = 0;
    let mut smiles_input = String::new();
    let mut smiles_popup_focus = true;
    let mut show_logs = false;
    let mut show_atom_popup = false;
    let mut show_neighbors = false;
    let mut show_ids = false;
    let collector = egui_tracing::EventCollector::new();
    Registry::default()
        .with(
            filter::Targets::new()
                .with_default(LevelFilter::INFO)
                .with_target("lute", LevelFilter::TRACE)
                .with_target("gui", LevelFilter::TRACE),
        )
        .with(collector.clone())
        .init();
    eframe::run_simple_native(
        "Lute GUI",
        eframe::NativeOptions::default(),
        move |ctx, _frame| {
            let focused = ctx.memory(|mem| mem.focused().is_some());
            let close = ctx.input_mut(|i| {
                let close = i.consume_shortcut(&egui::KeyboardShortcut::new(
                    egui::Modifiers::CTRL,
                    egui::Key::W,
                ));
                if i.consume_shortcut(&egui::KeyboardShortcut::new(
                    egui::Modifiers::CTRL,
                    egui::Key::D,
                )) {
                    show_logs = !show_logs;
                }
                if !focused {
                    if i.consume_key(Default::default(), egui::Key::A) {
                        show_atom_popup = !show_atom_popup;
                    }
                    if i.consume_key(Default::default(), egui::Key::N) {
                        show_neighbors = !show_neighbors;
                    }
                    if i.consume_key(Default::default(), egui::Key::I) {
                        show_ids = !show_ids;
                    }
                }
                close
            });
            if close {
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            }
            egui::SidePanel::left("sidebar")
                .resizable(true)
                .show(ctx, |ui| {
                    ui.horizontal_top(|ui| {
                        ui.label("Fragments");
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Min), |ui| {
                            if ui.button("C").clicked()
                                || ctx.input_mut(|i| {
                                    i.consume_shortcut(&egui::KeyboardShortcut::new(
                                        egui::Modifiers::CTRL,
                                        egui::Key::R,
                                    ))
                                })
                            {
                                arena = Arena::new();
                            }
                            {
                                let popup_id = ui.make_persistent_id("load_smiles");
                                let response = ui.button("+");
                                if response.clicked()
                                    || (ctx.memory(|mem| mem.focused().is_none())
                                        && ctx.input_mut(|i| {
                                            i.consume_key(Default::default(), egui::Key::L)
                                        }))
                                {
                                    ui.memory_mut(|mem| mem.open_popup(popup_id));
                                    smiles_popup_focus = true;
                                }
                                egui::popup_above_or_below_widget(
                                    ui,
                                    popup_id,
                                    &response,
                                    egui::AboveOrBelow::Below,
                                    egui::PopupCloseBehavior::IgnoreClicks,
                                    |ui| {
                                        ui.set_min_width(100.0);
                                        let editor = ui.add(
                                            egui::TextEdit::singleline(&mut smiles_input)
                                                .font(egui::TextStyle::Monospace)
                                                .return_key(None),
                                        );
                                        if smiles_popup_focus {
                                            editor.request_focus();
                                            smiles_popup_focus = false;
                                        }
                                        ui.with_layout(
                                            egui::Layout::top_down_justified(egui::Align::Center),
                                            |ui| {
                                                let load = ui.button("Load");
                                                let cancel = ui.button("Cancel");
                                                if load.clicked()
                                                    || (editor.has_focus()
                                                        && ctx.input_mut(|i| {
                                                            i.consume_key(
                                                                Default::default(),
                                                                egui::Key::Enter,
                                                            )
                                                        }))
                                                {
                                                    match SmilesParser::new(&smiles_input).parse() {
                                                        Ok(graph) => {
                                                            let idx = arena.insert_mol(&graph);
                                                            current = idx.0;
                                                        }
                                                        Err(err) => tracing::error!(
                                                            err = &err as &(dyn std::error::Error
                                                                  + Send
                                                                  + Sync
                                                                  + 'static),
                                                            "error parsing SMILES"
                                                        ),
                                                    }
                                                    smiles_input.clear();
                                                    ui.memory_mut(|mem| mem.close_popup());
                                                } else if cancel.clicked()
                                                    || (editor.has_focus()
                                                        && ctx.input_mut(|i| {
                                                            i.consume_key(
                                                                Default::default(),
                                                                egui::Key::Escape,
                                                            )
                                                        }))
                                                {
                                                    smiles_input.clear();
                                                    ui.memory_mut(|mem| mem.close_popup());
                                                }
                                            },
                                        );
                                    },
                                );
                            }
                        });
                    });
                    ui.separator();
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        ui.with_layout(egui::Layout::top_down_justified(egui::Align::Min), |ui| {
                            for n in 0..arena.expose_frags().len() {
                                let n = n as u32;
                                let button = ui.button(format!(
                                    "{n}: {}",
                                    arena.molecule(n.into()).smiles(Default::default())
                                ));
                                if button.clicked() {
                                    current = n;
                                }
                            }
                        });
                    });
                });
            if show_logs {
                egui::TopBottomPanel::bottom("logs")
                    .resizable(true)
                    .show(ctx, |ui| {
                        ui.add(egui_tracing::Logs::new(collector.clone()));
                    });
            }
            egui::CentralPanel::default().show(ctx, |ui| {
                if !arena.expose_frags().is_empty() {
                    let _guard =
                        tracing::subscriber::set_default(tracing::subscriber::NoSubscriber::new());
                    let mol = arena.molecule(current.into());
                    let click = ctx.pointer_interact_pos();
                    let mut selected = None;
                    ui.centered_and_justified(|ui| {
                        let mut atoms = mol
                            .node_references()
                            .map(|a| {
                                if a.weight().protons == 0 {
                                    85
                                } else {
                                    a.weight().protons
                                }
                            })
                            .collect::<Vec<_>>();
                        let mut edges = mol
                            .edge_references()
                            .map(|e| {
                                [
                                    e.source().0 as u16,
                                    e.target().0 as u16,
                                    e.weight().bond_count().floor() as u16,
                                ]
                            })
                            .collect::<Vec<_>>();
                        for atom in mol.node_references() {
                            let idx = atom.id().0 as u16;
                            let data = atom.weight().data;
                            for _ in 0..data.unknown() {
                                edges.push([idx as _, atoms.len() as _, 1]);
                                atoms.push(85);
                            }
                        }
                        let mut locs = coordgen::gen_coords(&atoms, &edges).unwrap();
                        let (min_x, min_y, max_x, max_y) = if atoms.is_empty() {
                            (0.0, 0.0, 0.0, 0.0)
                        } else {
                            locs.iter().fold(
                                (
                                    f32::INFINITY,
                                    f32::INFINITY,
                                    f32::NEG_INFINITY,
                                    f32::NEG_INFINITY,
                                ),
                                |(ix, iy, ax, ay), l| {
                                    (ix.min(l.0), iy.min(l.1), ax.max(l.0), ay.max(l.1))
                                },
                            )
                        };
                        let diff_x = max_x - min_x + 80.0;
                        let diff_y = max_y - min_y + 80.0;
                        let max_axis = diff_x.max(diff_y);

                        let painter = ui.painter();
                        let bounds = painter.clip_rect();
                        let size = (bounds.max - bounds.min).abs();
                        let min_axis = size.min_elem();
                        let scale = min_axis / max_axis * 0.75;
                        for (x, y) in &mut locs {
                            *x = (*x - min_x) * scale + bounds.min.x + (size.x - diff_x) / 2.0;
                            *y = (*y - min_y) * scale + bounds.min.y + (size.y - diff_y) / 2.0;
                        }

                        let mut neighbors = Vec::new();

                        for (aref, &(cx, cy)) in mol.node_references().zip(&locs) {
                            let atom = aref.weight();
                            let left =
                                [8, 9, 16, 17, 34, 35, 52, 53, 84, 85].contains(&atom.protons);
                            let is_protium = atom.protons == 1 && atom.isotope == 0;
                            let color = Color32::from_hex(svg_atom_color(atom.protons)).unwrap();
                            let rect = if let Some((n, (mut dx, mut dy))) = mol
                                .neighbors(aref.id())
                                .map(|i| (1usize, locs[i.0 as usize]))
                                .reduce(|(a, (ax, ay)), (b, (bx, by))| (a + b, (ax + bx, ay + by)))
                            {
                                dx /= n as f32;
                                dy /= n as f32;
                                dx -= cx;
                                dy -= cy;
                                let mag = (dx * dx + dy * dy).sqrt();
                                dx /= mag;
                                dy /= mag;
                                if show_ids {
                                    painter.text(
                                        egui::pos2(cx - dy * scale * 10.0, cy + dx * scale * 10.0),
                                        egui::Align2::CENTER_CENTER,
                                        aref.id().0,
                                        Default::default(),
                                        Color32::LIGHT_BLUE,
                                    );
                                }
                                if atom.protons != 6 || atom.isotope != 0 {
                                    let mut s = String::new();
                                    if dx > 0.0 && left {
                                        match atom.data.hydrogen() + is_protium as u8 {
                                            0 => {}
                                            1 => s.push('H'),
                                            n => {
                                                let _ = write!(s, "H{}", Subscript(n));
                                            }
                                        }
                                    }
                                    if atom.protons == 0 {
                                        s.push(match atom.isotope {
                                            0xFFFE => 'A',
                                            0xFFFC => 'Q',
                                            0x0100 => 'X',
                                            0x042B => 'M',
                                            _ => '*',
                                        });
                                    } else {
                                        if atom.isotope != 0 {
                                            let _ = write!(s, "{}", Superscript(atom.isotope));
                                        }
                                        if !is_protium {
                                            s.push_str(ATOM_DATA[atom.protons as usize].sym);
                                        }
                                    }
                                    if dx < 0.0 || !left {
                                        match atom.data.hydrogen() + is_protium as u8 {
                                            0 => {}
                                            1 => s.push('H'),
                                            n => {
                                                let _ = write!(s, "H{}", Subscript(n));
                                            }
                                        }
                                    }
                                    let mut rect = painter.text(
                                        egui::pos2(cx, cy),
                                        egui::Align2::CENTER_CENTER,
                                        s,
                                        FontId::proportional(scale * 10.0),
                                        color,
                                    );
                                    if atom.charge != 0 {
                                        let s = match atom.charge {
                                            1 => "+".to_string(),
                                            -1 => "-".to_string(),
                                            _ => format!("{:+}", atom.charge),
                                        };
                                        let sub = painter.text(
                                            egui::pos2(rect.max.x, rect.min.y),
                                            egui::Align2::LEFT_TOP,
                                            s,
                                            FontId::proportional(scale * 5.0),
                                            color,
                                        );
                                        rect = rect.union(sub);
                                    }
                                    rect
                                } else if atom.charge != 0 {
                                    painter
                                        .text(
                                            egui::pos2(cx - dx * scale, cy - dy * scale),
                                            egui::Align2::CENTER_CENTER,
                                            match atom.charge {
                                                1 => "+".to_string(),
                                                -1 => "-".to_string(),
                                                _ => format!("{:+}", atom.charge),
                                            },
                                            FontId::proportional(scale * 10.0),
                                            Color32::from_gray(0x60),
                                        )
                                        .union(egui::Rect::from_pos(egui::pos2(cx, cy)))
                                } else {
                                    egui::Rect::from_x_y_ranges(
                                        egui::Rangef::new(cx - scale * 5.0, cx + scale * 5.0),
                                        egui::Rangef::new(cy - scale * 5.0, cy + scale * 5.0),
                                    )
                                }
                            } else if atom.protons != 6 || atom.isotope != 0 {
                                let mut s = String::new();
                                if left {
                                    match atom.data.hydrogen() + is_protium as u8 {
                                        0 => {}
                                        1 => s.push('H'),
                                        n => {
                                            let _ = write!(s, "H{}", Subscript(n));
                                        }
                                    }
                                }
                                if atom.protons == 0 {
                                    s.push(match atom.isotope {
                                        0xFFFE => 'A',
                                        0xFFFC => 'Q',
                                        0x0100 => 'X',
                                        0x042B => 'M',
                                        _ => '*',
                                    });
                                } else {
                                    if atom.isotope != 0 {
                                        let _ = write!(s, "{}", Superscript(atom.isotope));
                                    }
                                    if !is_protium {
                                        s.push_str(ATOM_DATA[atom.protons as usize].sym);
                                    }
                                }
                                if !left {
                                    match atom.data.hydrogen() + is_protium as u8 {
                                        0 => {}
                                        1 => s.push('H'),
                                        n => {
                                            let _ = write!(s, "H{}", Subscript(n));
                                        }
                                    }
                                }
                                let mut rect = painter.text(
                                    egui::pos2(cx, cy),
                                    egui::Align2::CENTER_CENTER,
                                    s,
                                    FontId::proportional(scale * 10.0),
                                    color,
                                );
                                if atom.charge != 0 {
                                    let s = match atom.charge {
                                        1 => "+".to_string(),
                                        -1 => "-".to_string(),
                                        _ => format!("{:+}", atom.charge),
                                    };
                                    let sub = painter.text(
                                        egui::pos2(rect.max.x, rect.min.y),
                                        egui::Align2::LEFT_TOP,
                                        s,
                                        FontId::proportional(scale * 5.0),
                                        color,
                                    );
                                    rect = rect.union(sub);
                                }
                                rect
                            } else {
                                egui::Rect::from_x_y_ranges(
                                    egui::Rangef::new(cx - scale * 5.0, cx + scale * 5.0),
                                    egui::Rangef::new(cy - scale * 5.0, cy + scale * 5.0),
                                )
                            };
                            if let Some(pos) = click {
                                if rect.contains(pos) {
                                    selected = Some((aref.id(), rect));
                                    let center = rect.center();
                                    if show_atom_popup || show_neighbors {
                                        painter.circle_filled(
                                            center,
                                            rect.size().length() * 0.6,
                                            Color32::YELLOW.gamma_multiply(0.25),
                                        );
                                    }
                                    if show_neighbors {
                                        neighbors.extend(mol.edges(aref.id()));
                                    }
                                }
                            }
                        }

                        for edge in mol.edge_references() {
                            let [ix1, ix2] = std::cmp::minmax(
                                mol.to_index(edge.source()),
                                mol.to_index(edge.target()),
                            );
                            let (ox1, oy1) = locs[ix1];
                            let (ox2, oy2) = locs[ix2];
                            let (mut x1, mut y1) = (ox1, oy1);
                            let (mut x2, mut y2) = (ox2, oy2);
                            let (mut dx, mut dy) = (x2 - x1, y2 - y1);
                            let mag = (dx * dx + dy * dy).sqrt();
                            dx /= mag;
                            dy /= mag;
                            if atoms[ix1] != 6
                                || mol.node_weight(edge.source()).unwrap().isotope != 0
                            {
                                x1 += dx * scale * 9.0;
                                y1 += dy * scale * 9.0;
                            }
                            if atoms[ix2] != 6
                                || mol.node_weight(edge.target()).unwrap().isotope != 0
                            {
                                x2 -= dx * scale * 9.0;
                                y2 -= dy * scale * 9.0;
                            }
                            let (mut rdx, mut rdy) = (-dy * 2.75, dx * 2.75);
                            let stroke = if neighbors.iter().any(|e| e.id() == edge.id()) {
                                (4.0, Color32::GOLD)
                            } else {
                                (2.0, Color32::LIGHT_GRAY)
                            };
                            match *edge.weight() {
                                Bond::Non => {}
                                Bond::Single | Bond::Left | Bond::Right => {
                                    painter.line_segment(
                                        [egui::pos2(x1, y1), egui::pos2(x2, y2)],
                                        stroke,
                                    );
                                }
                                Bond::Double | Bond::DoubleE | Bond::DoubleZ => {
                                    painter.line_segment(
                                        [
                                            egui::pos2(x1 - rdx, y1 - rdy),
                                            egui::pos2(x2 - rdx, y2 - rdy),
                                        ],
                                        stroke,
                                    );
                                    painter.line_segment(
                                        [
                                            egui::pos2(x1 + rdx, y1 + rdy),
                                            egui::pos2(x2 + rdx, y2 + rdy),
                                        ],
                                        stroke,
                                    );
                                }
                                Bond::Triple => {
                                    painter.line_segment(
                                        [
                                            egui::pos2(x1 - 2.0 * rdx, y1 - 2.0 * rdy),
                                            egui::pos2(x2 - 2.0 * rdx, y2 - 2.0 * rdy),
                                        ],
                                        stroke,
                                    );
                                    painter.line_segment(
                                        [
                                            egui::pos2(x1 + 2.0 * rdx, y1 + 2.0 * rdy),
                                            egui::pos2(x2 + 2.0 * rdx, y2 + 2.0 * rdy),
                                        ],
                                        stroke,
                                    );
                                    painter.line_segment(
                                        [egui::pos2(x1, y1), egui::pos2(x2, y2)],
                                        stroke,
                                    );
                                }
                                Bond::Quad => {
                                    painter.line_segment(
                                        [
                                            egui::pos2(x1 - 3.0 * rdx, y1 - 3.0 * rdy),
                                            egui::pos2(x2 - 3.0 * rdx, y2 - 3.0 * rdy),
                                        ],
                                        stroke,
                                    );
                                    painter.line_segment(
                                        [
                                            egui::pos2(x1 - rdx, y1 - rdy),
                                            egui::pos2(x2 - rdx, y2 - rdy),
                                        ],
                                        stroke,
                                    );
                                    painter.line_segment(
                                        [
                                            egui::pos2(x1 + rdx, y1 + rdy),
                                            egui::pos2(x2 + rdx, y2 + rdy),
                                        ],
                                        stroke,
                                    );
                                    painter.line_segment(
                                        [
                                            egui::pos2(x1 + 3.0 * rdx, y1 + 3.0 * rdy),
                                            egui::pos2(x2 + 3.0 * rdx, y2 + 3.0 * rdy),
                                        ],
                                        stroke,
                                    );
                                }
                                Bond::Aromatic => {
                                    let mut flip = 0.0;
                                    let id1 = edge.source();
                                    let id2 = edge.target();
                                    for e in mol.edges(id1).chain(mol.edges(id2)) {
                                        if *e.weight() != Bond::Aromatic {
                                            continue;
                                        }
                                        if [id1, id2].contains(&e.target()) {
                                            continue;
                                        }
                                        let (x1, y1) = locs[e.source().0 as usize];
                                        let (x2, y2) = locs[e.target().0 as usize];
                                        let ndx = x2 - x1;
                                        let ndy = y2 - y1;
                                        let dot = dx * ndy - dy * ndx;
                                        flip += dot;
                                    }
                                    if flip > -0.0001 {
                                        rdx = -rdx;
                                        rdy = -rdy;
                                    }
                                    rdx *= 1.5;
                                    rdy *= 1.5;
                                    painter.add(egui::Shape::dashed_line(
                                        &[
                                            egui::pos2(x1 - rdx, y1 - rdy),
                                            egui::pos2(x2 - rdx, y2 - rdy),
                                        ] as _,
                                        stroke,
                                        10.0,
                                        10.0,
                                    ));
                                    painter.line_segment(
                                        [egui::pos2(x1, y1), egui::pos2(x2, y2)],
                                        stroke,
                                    );
                                }
                                _ => panic!("invalid bond!"),
                            }
                        }

                        let mut idx = mol.node_count();

                        for (atom, &(x1, y1)) in mol.node_references().zip(&locs) {
                            let atom = atom.weight();
                            let data = atom.data;
                            for i in 0..data.unknown() {
                                let (x2, y2) = locs[idx + (i as usize)];
                                let mut dx = x2 - x1;
                                let mut dy = y2 - y1;
                                let mul = (dx * dx + dy * dy).sqrt().recip() * scale * 12.0;
                                dx *= mul;
                                dy *= mul;
                                let x3 = x2 - dx;
                                let y3 = y2 - dy;
                                let mut x0 = x1;
                                let mut y0 = y1;
                                if atom.protons != 6 || atom.isotope != 0 {
                                    x0 += dx;
                                    y0 += dy;
                                }
                                painter.line_segment(
                                    [egui::pos2(x0, y0), egui::pos2(x3, y3)],
                                    (2.0, Color32::LIGHT_GRAY),
                                );
                                painter.text(
                                    egui::pos2(x2, y2),
                                    egui::Align2::CENTER_CENTER,
                                    "R",
                                    FontId::proportional(scale * 15.0),
                                    Color32::from_rgb(0x40, 0x7F, 0x00),
                                );
                            }
                            idx += data.unknown() as usize;
                        }
                    });

                    if show_atom_popup {
                        if let Some((id, rect)) = selected {
                            ui.allocate_new_ui(
                                egui::UiBuilder::new().max_rect(egui::Rect::from_min_size(
                                    rect.center(),
                                    egui::vec2(100.0, 100.0),
                                )),
                                |ui| {
                                    egui::Frame::popup(ui.style()).show(ui, |ui| {
                                        let atom = mol.node_weight(id).unwrap();
                                        ui.label(format!("{:#}", atom));
                                        ui.label(format!("protons: {}", atom.isotope));
                                        ui.label(format!("isotope: {}", atom.isotope));
                                        ui.separator();
                                        ui.label(format!("hydrogen: {}", atom.data.hydrogen()));
                                        ui.label(format!("single bonds: {}", atom.data.single()));
                                        ui.label(format!("placeholders: {}", atom.data.unknown()));
                                        ui.label(format!("other bonds: {}", atom.data.other()));
                                    });
                                },
                            );
                        }
                    }

                    ui.with_layout(egui::Layout::bottom_up(egui::Align::Max), |ui| {
                        egui::Frame::group(ui.style()).show(ui, |ui| {
                            let plain = egui::TextFormat::default();
                            let code = egui::TextFormat {
                                background: ui.style().visuals.code_bg_color,
                                ..egui::TextFormat::simple(
                                    FontId::monospace(14.0),
                                    Color32::PLACEHOLDER,
                                )
                            };
                            let mut job = egui::text::LayoutJob::single_section(
                                format!("Atomic mass: {:.3}u\nFast SMILES: ", mol.mass()),
                                plain.clone(),
                            );
                            job.append(
                                &mol.smiles(SmilesConfig::fast_roundtrip()),
                                4.0,
                                code.clone(),
                            );
                            job.append("\nCanon SMILES: ", 0.0, plain);
                            job.append(&mol.smiles(SmilesConfig::new()), 4.0, code);
                            ui.label(job);
                        });
                    });
                }
            });
        },
    )
}
